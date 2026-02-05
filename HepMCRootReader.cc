

#include "HepMCRootReader.hh"

#include "GenEventData.hh"

#include "TFile.h"
#include "TTree.h"

#include <iostream>


HepMCRootReader::HepMCRootReader( const std::string & filename,
				  const std::string & treename,
				  const std::string & branchname,
				  const Double_t mcnonradcutin ) :
  nevents(0), hepmcrootfile(nullptr), hepmctree(nullptr), mcnonradcut(mcnonradcutin),
  eventData(nullptr), runInfoData(nullptr),
  vtlvCache{ { "parton", std::vector<TLorentzVector>() },
    { "hadron", std::vector<TLorentzVector>() } },
  cacheIsValid{ { "parton", false }, { "hadron", false } } {
  std::cout << "HepMCRootReader::HepMCRootReader: opening " << filename << std::endl;
  hepmcrootfile= TFile::Open( filename.c_str() );
  if( not hepmcrootfile ) {
    std::string txt= "HepMCRootReader::HepMCRootReader: file not found: ";
    txt+= filename;
    throw std::runtime_error( txt );
  }
  if( not hepmcrootfile->IsOpen() ) {
    std::string txt= "HepMCRootReader::HepMCRootReader: file not open: ";
    txt+= filename;
    throw std::runtime_error( txt );
  }
  hepmctree= reinterpret_cast<TTree*>( hepmcrootfile->Get( treename.c_str() ) );
  if( hepmctree == nullptr ) {
    std::string txt= "HepMCRootReader::HepMCRootReader: tree not found: ";
    txt+= treename;
    throw std::runtime_error( txt );
  }
  std::cout << "HepMCRootReader::HepMCRootReader: " 
	    << GetNumberEntries() << " events on file" << std::endl;
  eventData= new HepMC3::GenEventData();
  Int_t status= hepmctree->SetBranchAddress( branchname.c_str(), &eventData );
  if( status < 0 ) {
    throw std::runtime_error( "HepMCRootReader::HepMCRootReader: can't create event branch hepmc3_event" );
  }
  runInfoData= reinterpret_cast<HepMC3::GenRunInfoData*>( hepmcrootfile->Get( "GenRunInfoData" ) );
  if( runInfoData == nullptr ) {
    std::cout << "HepMCRootReader::HepMCRootReader: no GenRunInfoData" << std::endl;
  }
  
}

HepMCRootReader::~HepMCRootReader() {
  std::cout << "HepMCRootReader::~HepMCRootReader: closing file "
	    << hepmcrootfile->GetName() << std::endl;
  hepmcrootfile->Close();
  delete hepmcrootfile;
  delete eventData;
}

Int_t HepMCRootReader::GetNumberEntries() {
  if( hepmctree ) return hepmctree->GetEntries();
  else throw std::runtime_error( "HepMCRootReader::GetNumberEntries: nullptr" );
}

// Helpers to navigate event record
// Compare object addresses for object quality
template <typename T>
class CompareObjectAddress {
  const T* ref;
public:
  CompareObjectAddress( const T* obj ) : ref(obj) {}
  bool operator()( const T & obj ) {
    return &obj == ref;
  }
};
// Find position (index) in vector
template <typename T>
size_t findInVector( const std::vector<T> & v, const T* obj ) {
  CompareObjectAddress<T> coa( obj );
  auto it= std::find_if( v.begin(), v.end(), coa );
  auto index= std::distance( v.begin(), it );
  return index;
}

bool HepMCRootReader::GetEvent( Int_t ievent ) {  
  bool result= false;
  if( hepmctree and hepmctree->GetEvent( ievent ) > 0 ) {
    if( ievent % 100 == 0 ) {
	std::cout << "HepMCRootReader::GetEvent: event " << ievent << std::endl;
    }
    result= true;
    for( const auto & keyValue : cacheIsValid ) cacheIsValid[keyValue.first]= false;
    nevents++;
    incomingParticlesMap.clear();
    outgoingParticlesMap.clear();
    beginVertexMap.clear();
    endVertexMap.clear();
    // Setup links between particles and vertices in maps:
    for( size_t i= 0; i < eventData->links1.size(); i++ ) {
      Int_t id1= eventData->links1[i];
      Int_t id2= eventData->links2[i];
      if( id1 > 0 ) {
	const Vertex & vertex= eventData->vertices[ -id2-1 ];
	ParticleVector & incoming= incomingParticlesMap[ &vertex ];
	const Particle & particle= eventData->particles[ id1-1 ];
	incoming.push_back( &particle );
	endVertexMap[ &particle ]= &vertex;
      }
      else {
	const Vertex & vertex= eventData->vertices[ -id1-1 ];
	ParticleVector & outgoing= outgoingParticlesMap[ &vertex ];
	const Particle & particle= eventData->particles[ id2-1 ];
	outgoing.push_back( &particle );
	beginVertexMap[ &particle ]= &vertex;
      }
    }
    // ISR photons
    findISRphotons();
  }
  else {
    std::cout << "HepMCRootReader::GetEvent: event " << nevents << " failure, stop"
	      << std::endl;
    result= false;
  }
  return result;
}
bool HepMCRootReader::GetNextEvent( Int_t maxevt ) {
  bool maxreached= maxevt > 0 and maxevt == nevents;
  return not maxreached and GetEvent( nevents );
}

// Get connected particles from a vertex
HepMCRootReader::ParticleVector
getParticles( const HepMCRootReader::Particle* particle,
	      const HepMCRootReader::ParticleVertexMap & pvm,
	      const HepMCRootReader::VertexParticlesMap & vpm ) {
  if( not ( pvm.count( particle ) == 0 ) ) {
    const HepMCRootReader::Vertex* vertex= pvm.at( particle );
    if( vertex ) {
      const HepMCRootReader::ParticleVector& particles= vpm.at( vertex );
      return particles;
    }
  }
  return HepMCRootReader::ParticleVector();
}
// Selector interface for event iterator
class ParticleSelector {
public:
  virtual bool operator()( const HepMCRootReader::Particle* ) const { return true; }
};
// Event iterator to move down or up event history
class EventIterator {
  const HepMCRootReader::Particle* particle;
  const HepMCRootReader::ParticleVertexMap vertexMap;
  const HepMCRootReader::VertexParticlesMap particlesMap;
  const ParticleSelector & selector;
public:
  EventIterator( const HepMCRootReader::Particle* p,
		 const HepMCRootReader::ParticleVertexMap & pvm,
		 const HepMCRootReader::VertexParticlesMap & vpm,
		 const ParticleSelector & ps=ParticleSelector() ) :
    particle(p), vertexMap(pvm), particlesMap(vpm), selector(ps) {}
  bool next() {
    const HepMCRootReader::ParticleVector& vertexParticles=
      getParticles( particle, vertexMap, particlesMap );
    for( const HepMCRootReader::Particle* vertexParticle : vertexParticles ) {
      if( selector( vertexParticle ) ) {
	particle= vertexParticle;
	return true;
      }
    }
    return false;
  }
  const HepMCRootReader::Particle* getParticle() { return particle; }
  const HepMCRootReader::Vertex* getVertex() { return vertexMap.at( particle ); }
};
// Select identical particles with same PID
class IdenticalParticleSelector : public ParticleSelector {
  const HepMCRootReader::Particle* particle;
public:
  IdenticalParticleSelector( const HepMCRootReader::Particle* p ) :
    particle(p) {}
  virtual bool operator()( const HepMCRootReader::Particle* p ) const {    
    return ( particle->pid == p->pid and p->status != 3 );
  }
};

// Photons from initial beam particles
void HepMCRootReader::findISRphotons() {
  ISRphotons.clear();
  for( const Particle* beamParticle : getBeamParticles() ) {
    // for each beam particle go down history and find photons
    IdenticalParticleSelector ips( beamParticle );
    EventIterator eventIter( beamParticle, endVertexMap, outgoingParticlesMap,
			     ips );
    while( eventIter.next() ) {
      const ParticleVector outgoing= getDaughters( eventIter.getParticle() );
      for( const Particle* outparticle : outgoing ) {
	if( outparticle->pid == 22 and outparticle->status == 1 ) {
	  ISRphotons.push_back( outparticle );
	}
      }
    }
  }
  return;
}

HepMCRootReader::ParticleVector HepMCRootReader::getBeamParticles() {
  // Initial beam particles "status 4"
  ParticleVector beamParticles;
  for( const Particle & particle : eventData->particles ) {
    if( particle.status == 4 ) beamParticles.push_back( &particle );
  }
  if( beamParticles.size() != 2 ) {
    throw std::runtime_error( "HepMCRootReader::getBeamParticles: no 2 beam particles" );
  }
  return beamParticles;
}

Double_t HepMCRootReader::event_scale() {
  TLorentzVector sum;
  for( const Particle* beamParticle : getBeamParticles() ) {
    TLorentzVector tlv( beamParticle->momentum.m_v1,
			beamParticle->momentum.m_v2,
			beamParticle->momentum.m_v3,
			beamParticle->momentum.m_v4 );
    sum+= tlv;
  }
  return sum.M();
}  

const std::vector<TLorentzVector>
HepMCRootReader::GetLorentzVectors( const std::string & opt ) {
  if( cacheIsValid[opt] ) return vtlvCache[opt];
  if( opt == "hadron" ) {
    getHadron();
  }
  else if( opt == "parton" ) {
    getParton();
  }
  else if( opt == "isr" ) {
    getIsr();
  }
  vtlvCache[opt]= vtlv;
  cacheIsValid[opt]= true;
  return vtlv;
}

void HepMCRootReader::getIsr() {
  vtlv.clear();
  // ISRphotons filled in getEvent()
  for( const Particle* photon : ISRphotons ) {
    TLorentzVector tlv( photon->momentum.m_v1,
			photon->momentum.m_v2,
			photon->momentum.m_v3,
			photon->momentum.m_v4 );
    vtlv.push_back( tlv );
  }
  return;
}

bool isParton( const HepMCRootReader::Particle* particle ) {
  // Quarks, leptons, photon, gluon, W, Z
  bool result= false;
  Int_t abspid= abs( particle->pid );
  if( abspid <= 6 or
      ( abspid >= 11 and abspid <= 16 ) or
      ( abspid >= 21 and abspid <= 24 ) ) {
    result= true;
  }
  return result;
}

void HepMCRootReader::getParton() {
  vtlv.clear();
  TLorentzVector sum;
  Int_t ip= -1;
  std::vector<Int_t> selected;

  //std::cout << "HepMCRootReader::getParton: print event record" << std::endl;
  //printParticlesVertices();
  //std::cout << "HepMCRootReader::getParton: before particles loop" << std::endl;
  
  for( const Particle & particle : eventData->particles ) {
    ip++;

    //std::cout << "HepMCRootReader::getParton: in particles loop ip= " << ip << std::endl;
    
    if( not isParton( &particle ) ) {
      //std::cout << "HepMCRootReader::getParton: not a parton ip= " << ip << std::endl;
      continue;
    }
    // Detect decays of partons to partons:
    if( endVertexMap.count( &particle ) == 0 ) {
      //std::cout << "HepMCRootReader::getParton: not in endVertexMap ip= " << ip << std::endl;
      continue;
    }
    // const Vertex* endVertex= endVertexMap.at( &particle );
    // const ParticleVector & decayProducts= outgoingParticles.at( endVertex );
    ParticleVector decayProducts= getDaughters( &particle );
    bool partonDecay= true;
    for( const Particle* decayParticle : decayProducts ) {
      if( not isParton( decayParticle ) ) {
	partonDecay= false;
	break;
      }
    }
    if( partonDecay ) {
      //std::cout << "HepMCRootReader::getParton: decay to partons ip= " << ip << std::endl;
      continue;
    }
    // Check the parton does not come from a partonic hadron decay
    ParticleVector parents= getParents( &particle );
    bool partonicHadronDecay= false;
    for( const Particle* parent : parents ) {
      if( abs( parent->pid ) > 100 ) {
	partonicHadronDecay= true;
	break;
      }
    }
    if( partonicHadronDecay ) {
      //std::cout << "HepMCRootReader::getParton: from partonic hadron decay ip= " << ip
      //	<< std::endl;
      continue;
    }
    // Selected parton
    selected.push_back( ip );
    TLorentzVector tlv( particle.momentum.m_v1,
			particle.momentum.m_v2,
			particle.momentum.m_v3,
			particle.momentum.m_v4 );
    vtlv.push_back( tlv );
    sum+= tlv;
  }
  // Check 4-momentum conservation
  for( const Particle* isrphoton : ISRphotons ) {
    TLorentzVector tlv( isrphoton->momentum.m_v1,
			isrphoton->momentum.m_v2,
			isrphoton->momentum.m_v3,
			isrphoton->momentum.m_v4 );
    sum+= tlv;
  }
  if( fabs( sum.Px() ) > 1.0e-6 or
      fabs( sum.Py() ) > 1.0e-6 or
      fabs( sum.Pz() ) > 1.0e-6 or
      fabs( sum.E() - event_scale() ) > 1.0e-6 ) {
    std::cout << "HepMCRootReader::getParton: 4-vector sum mismatch" << std::endl;
    std::cout << "4-momentum sum: "
	      << sum.Px() << " "
	      << sum.Py() << " "
	      << sum.Pz() << " "
	      << sum.E() << ", event scale: "
	      << event_scale() << std::endl;
    std::cout << "Selected: ";
    for( Int_t ip : selected ) std::cout << ip << " ";
    std::cout << std::endl;
    printParticlesVertices();
  }
  return;
}

void HepMCRootReader::getHadron() {
  vtlv.clear();
  TLorentzVector sum;
  for( const Particle & particle : eventData->particles ) {
    // Stable particle and not tagged as ISR
    if( particle.status == 1 and
	std::find( ISRphotons.begin(), ISRphotons.end(),
		   &particle ) == ISRphotons.end() ) {
      TLorentzVector tlv( particle.momentum.m_v1,
			  particle.momentum.m_v2,
			  particle.momentum.m_v3,
			  particle.momentum.m_v4 );
      vtlv.push_back( tlv );
      sum+= tlv;
    }
  }
  // Check 4-momentum conservation
  for( const Particle* isrphoton : ISRphotons ) {
    TLorentzVector tlv( isrphoton->momentum.m_v1,
			isrphoton->momentum.m_v2,
			isrphoton->momentum.m_v3,
			isrphoton->momentum.m_v4 );   
    sum+= tlv;
  }
  if( fabs( sum.Px() ) > 1.0e-6 or
      fabs( sum.Py() ) > 1.0e-6 or
      fabs( sum.Pz() ) > 1.0e-6 or
      fabs( sum.E() - event_scale() ) > 1.0e-6 ) {
    std::cout << "HepMCRootReader::getHadron: 4-vector sum mismatch event "
	      << eventData->event_number
	      << std::endl;
    std::cout << "4-momentum sum: "
	      << sum.Px() << " "
	      << sum.Py() << " "
	      << sum.Pz() << " "
	      << sum.E() << " event scale: "
	      << event_scale() << ", difference: "
	      << sum.E() - event_scale()
	      << std::endl;
    printParticlesVertices();
  }
  return;
}

const std::map<std::string,bool>
HepMCRootReader::getSelections( const std::string & opt ) {
  std::map<std::string,bool> result;
  return result;
}

bool HepMCRootReader::MCNonRad() {
  Double_t sumE= 0.0;
  for( const Particle* isrphoton : ISRphotons ) {
    sumE+= isrphoton->momentum.m_v4;
  }
  bool result= true;
  if( sumE > mcnonradcut ) result= false;
  cutflow["mcnonrad"]= result;
  return result;
}

bool HepMCRootReader::isMC() {
  return true;
}

void HepMCRootReader::printParticlesVertices() {
  Int_t ip= 0;
  std::cout << "HepMCRootReader::printParticlesVertices: GenEvent particles" << std::endl;
  for( const Particle & particle : eventData->particles ) {
    std::cout << ip << " " << particle.pid << " " << particle.status << " "
	      << particle.mass << " "
	      << particle.momentum.m_v1 << " "
	      << particle.momentum.m_v2 << " "
	      << particle.momentum.m_v3 << " "
	      << particle.momentum.m_v4 << " ";
    std::cout << "parents ";
    ParticleVector parents= getParents( &particle );
    std::vector<size_t> parentIDs= getIndices( parents );
    for( const size_t index : parentIDs ) std::cout << index << " ";
    std::cout << "daughters ";
    ParticleVector daughters= getDaughters(  &particle );
    std::vector<size_t> daughterIDs= getIndices( daughters );
    for( const size_t index : daughterIDs ) std::cout << index << " ";
    std::cout << std::endl;
    ip++;
  }
  std::cout << "HepMCRootReader::printParticlesVertices: GenEvent ISR photons" << std::endl;
  for( const Particle* particle : ISRphotons ) {
    size_t index= findInVector<Particle>( eventData->particles, particle );
    std::cout << index << " " << particle->pid << " " << particle->status << " "
	      << particle->mass << " "
	      << particle->momentum.m_v1 << " "
	      << particle->momentum.m_v2 << " "
	      << particle->momentum.m_v3 << " "
	      << particle->momentum.m_v4 << " ";
    std::cout << "parents ";
    ParticleVector parents= getParents( particle );
    std::vector<size_t> parentIDs= getIndices( parents );
    for( const size_t index : parentIDs ) std::cout << index << " ";
    std::cout << "daughters ";
    ParticleVector daughters= getDaughters( particle );
    std::vector<size_t> daughterIDs= getIndices( daughters );
    for( const size_t index : daughterIDs ) std::cout << index << " ";
    std::cout << std::endl;
  }
  std::cout << "HepMCRootReader::printParticlesVertices: GenEvent vertices" << std::endl;
  for( size_t iv= 0; iv < eventData->vertices.size(); iv++ ) {
    const Vertex* vertex= &(eventData->vertices[iv]);
    std::cout << "vertex " << iv;
    if( incomingParticlesMap.count( vertex ) != 0 ) {
      std::cout << " incoming";
      const ParticleVector & incoming= incomingParticlesMap.at( vertex );
      for( const Particle* particle : incoming ) {
	size_t index= findInVector<Particle>( eventData->particles, particle );
	if( index == eventData->particles.size() ) {
	  std::cout << " particle not in eventData->particles ";
	}
	else {
	  std::cout << " " << index;
	}
      }
    }
    else {
      std::cout << " not in incoming";
    }
    if( outgoingParticlesMap.count( vertex ) != 0 ) {
      std::cout << " outgoing ";
      const ParticleVector & outgoing= outgoingParticlesMap.at( vertex );
      for( const Particle* particle : outgoing ) {
	size_t index= findInVector<Particle>( eventData->particles, particle );
	if( index == eventData->particles.size() ) {
	  std::cout << " particle not in eventData->particles ";
	}
	else {
	  std::cout << " " << index;
	}
      }      
    }
    std::cout << std::endl;
  }
  return;
}

// Particle parents and daughters
std::vector<size_t>
HepMCRootReader::getIndices( const ParticleVector& particles ) const {
  std::vector<size_t> result;
  for( const Particle* particle : particles ) {
    size_t index= findInVector<Particle>( eventData->particles, particle );
    result.push_back( index );
  }
  return result;
}
HepMCRootReader::ParticleVector
HepMCRootReader::getParents( const HepMCRootReader::Particle* particle ) const {
  const HepMCRootReader::ParticleVector& parents= getParticles( particle,
								beginVertexMap,
								incomingParticlesMap );
  return parents;
}
HepMCRootReader::ParticleVector
HepMCRootReader::getDaughters( const HepMCRootReader::Particle* particle ) const {
  const HepMCRootReader::ParticleVector& daughters= getParticles( particle,
								  endVertexMap,
								  outgoingParticlesMap );
  return daughters;
}


