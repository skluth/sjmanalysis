

#include "HepMCRootReader.hh"

#include "GenEventData.hh"

#include "TFile.h"
#include "TTree.h"

#include <iostream>


HepMCRootReader::HepMCRootReader( const std::string & filename,
				  const std::string & treename,
				  const std::string & branchname ) :
  nevents(0), hepmcrootfile(nullptr), hepmctree(nullptr),
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
  eventData= new HepMC::GenEventData();
  Int_t status= hepmctree->SetBranchAddress( branchname.c_str(), &eventData );
  if( status < 0 ) {
    throw std::runtime_error( "HepMCRootReader::HepMCRootReader: can't create event branch hepmc3_event" );
  }
  runInfoData= reinterpret_cast<HepMC::GenRunInfoData*>( hepmcrootfile->Get( "GenRunInfoData" ) );
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

// Helpers to navigate event record: compare object addresses for object quality
template <typename T>
class CompareObjectAddress {
  const T* ref;
public:
  CompareObjectAddress( const T* obj ) : ref(obj) {}
  bool operator()( const T & obj ) {
    return &obj == ref;
  }
};
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
    result= true;
    for( const auto & keyValue : cacheIsValid ) cacheIsValid[keyValue.first]= false;
    nevents++;
    incomingParticles.clear();
    outgoingParticles.clear();
    beginVertexMap.clear();
    endVertexMap.clear();
    // Setup links between particles and vertices in maps:
    for( size_t i= 0; i < eventData->links1.size(); i++ ) {
      Int_t id1= eventData->links1[i];
      Int_t id2= eventData->links2[i];
      if( id1 > 0 ) {
	const Vertex & vertex= eventData->vertices[ -id2-1 ];
	ParticleVector & incoming= incomingParticles[ &vertex ];
	const Particle & particle= eventData->particles[ id1-1 ];
	incoming.push_back( &particle );
	endVertexMap[ &particle ]= &vertex;
      }
      else {
	const Vertex & vertex= eventData->vertices[ -id1-1 ];
	ParticleVector & outgoing= outgoingParticles[ &vertex ];
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

class ParticleSelector {
public:
  virtual bool operator()( const HepMCRootReader::Particle* ) const { return true; }
};
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
    if( vertexMap.count( particle ) == 0 ) return false;
    const HepMCRootReader::Vertex* vertex= vertexMap.at( particle );
    if( particlesMap.count( vertex ) == 0 ) return false;
    const HepMCRootReader::ParticleVector & vertexParticles= particlesMap.at( vertex );
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

class IdenticalParticleSelector : public ParticleSelector {
  const HepMCRootReader::Particle* particle;
public:
  IdenticalParticleSelector( const HepMCRootReader::Particle* p ) :
    particle(p) {}
  virtual bool operator()( const HepMCRootReader::Particle* p ) const {    
    return ( particle->pid == p->pid and p->status != 3 );
  }
};
void HepMCRootReader::findISRphotons() {
  ISRphotons.clear();
  for( const Particle* beamParticle : getBeamParticles() ) {
    IdenticalParticleSelector ips( beamParticle );
    EventIterator eventIter( beamParticle, endVertexMap, outgoingParticles,
			     ips );
    while( eventIter.next() ) {      
      const Vertex* vertex= eventIter.getVertex();
      const ParticleVector & outgoing= outgoingParticles.at( vertex );
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


  printParticlesVertices();
  
  
  for( const Particle & particle : eventData->particles ) {
    ip++;
    if( not isParton( &particle ) ) continue;
    // Detect decays to partons:
    if( endVertexMap.count( &particle ) == 0 ) continue;
    const Vertex* endVertex= endVertexMap.at( &particle );
    const ParticleVector & decayProducts= outgoingParticles.at( endVertex );
    bool partonDecay= true;
    for( const Particle* decayParticle : decayProducts ) {
      if( not isParton( decayParticle ) ) {
	partonDecay= false;
	break;
      }
    }
    if( partonDecay ) continue;
    // Check the parton does not come from a partonic hadron decay
    bool partonicHadronDecay= false;
    EventIterator eventIter( &particle,
			     beginVertexMap,
			     incomingParticles );
    while( eventIter.next() ) {
      const HepMCRootReader::Particle* nextParticle= eventIter.getParticle();
      if( abs( nextParticle->pid ) > 100 ) {
	partonicHadronDecay= true;
	break;
      }
    }    
    if( partonicHadronDecay ) continue;
    selected.push_back( ip );
    TLorentzVector tlv( particle.momentum.m_v1,
			particle.momentum.m_v2,
			particle.momentum.m_v3,
			particle.momentum.m_v4 );
    vtlv.push_back( tlv );
    sum+= tlv;
  }
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
    printParticlesVertices();
    std::cout << "4-momentum sum: "
	      << sum.Px() << " "
	      << sum.Py() << " "
	      << sum.Pz() << " "
	      << sum.E() << ", event scale: "
	      << event_scale() << std::endl;
    std::cout << "Selected: ";
    for( Int_t ip : selected ) std::cout << ip << " ";
    std::cout << std::endl;
    std::cout << "HepMCRootReader::getParton2: 4-vector sum mismatch" << std::endl;
  }
  return;
}

void HepMCRootReader::getHadron() {
  vtlv.clear();
  TLorentzVector sum;
  for( const Particle & particle : eventData->particles ) {
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
    printParticlesVertices();
    std::cout << "4-momentum sum: "
	      << sum.Px() << " "
	      << sum.Py() << " "
	      << sum.Pz() << " "
	      << sum.E() << " event scale: "
	      << event_scale() << ", difference: "
	      << sum.E() - event_scale()
	      << std::endl;
    std::cout << "HepMCRootReader::getHadron: 4-vector sum mismatch event "
	      << eventData->event_number
	      << std::endl;
  }
  return;
}

const std::map<std::string,bool>
HepMCRootReader::getSelections( const std::string & opt ) {
  std::map<std::string,bool> result;
  return result;
}

bool HepMCRootReader::MCNonRad() {
  return true;
}


bool HepMCRootReader::isMC() {
  return true;
}

void HepMCRootReader::printParticlesVertices() {
  Int_t ip= 0;
  for( const Particle & particle : eventData->particles ) {
    std::cout << ip << " " << particle.pid << " " << particle.status << " "
	      << particle.mass << " "
	      << particle.momentum.m_v1 << " "
	      << particle.momentum.m_v2 << " "
	      << particle.momentum.m_v3 << " "
	      << particle.momentum.m_v4 << " "
	      << std::endl;
    ip++;
  }
  std::cout << "ISR photons" << std::endl;
  for( const Particle* particle : ISRphotons ) {
    size_t index= findInVector<Particle>( eventData->particles, particle );
    std::cout << index << " " << particle->pid << " " << particle->status << " "
	      << particle->mass << " "
	      << particle->momentum.m_v1 << " "
	      << particle->momentum.m_v2 << " "
	      << particle->momentum.m_v3 << " "
	      << particle->momentum.m_v4 << " "
	      << std::endl;
  }
  for( size_t iv= 0; iv < eventData->vertices.size(); iv++ ) {
    const Vertex* vertex= &(eventData->vertices[iv]);
    std::cout << "vertex " << iv;
    if( incomingParticles.count( vertex ) != 0 ) {
      std::cout << " incoming";
      const ParticleVector & incoming= incomingParticles.at( vertex );
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
    if( outgoingParticles.count( vertex ) != 0 ) {
      std::cout << " outgoing ";
      const ParticleVector & outgoing= outgoingParticles.at( vertex );
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

