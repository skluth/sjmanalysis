
#include "HepMC2Reader.hh"

#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"
#include "HepMC/IO_Exception.h"

#include <iostream>
#include <vector>
#include <algorithm>

HepMC2Reader::HepMC2Reader( const std::string & filename ) :
  input( filename ), nevent(0),
  vtlvCache{ { "parton", std::vector<TLorentzVector>() },
    { "hadron", std::vector<TLorentzVector>() } },
  cacheIsValid{ { "parton", false }, { "hadron", false } } {
  std::cout << "HepMC2Reader::HepMC2Reader: opened hepmc2 file " << filename << std::endl;
}

HepMC2Reader::~HepMC2Reader() {
  std::cout << "HepMC2Reader::~HepMC2Reader: processed events: " << nevent << std::endl;
}

// Follow a particle down until its end:
class ParticleHistory {
  const HepMC::GenParticle* particle;
public:
  ParticleHistory( const HepMC::GenParticle* p ) : particle(p) {}
  bool nextVertex() {
    bool result= false;
    const HepMC::GenVertex* vertex= particle->end_vertex();
    for( HepMC::GenVertex::particles_out_const_iterator
	   poutit= vertex->particles_out_const_begin();
	 poutit != vertex->particles_out_const_end(); poutit++ ) {
      const HepMC::GenParticle* outparticle= *poutit;
      if( particle->pdg_id() == outparticle->pdg_id() and outparticle->status() != 3 ) {
	particle= outparticle;
	result= true;
	break;
      }
    }
    return result;
  }
  const HepMC::GenVertex* getVertex() { return particle->end_vertex(); }
};
// Find ISR from beams:
void HepMC2Reader::findISRphotons() {
  std::pair<HepMC::GenParticle*,HepMC::GenParticle*> incoming= event.beam_particles();
  ISRphotons.clear();
  for( const HepMC::GenParticle* beamParticle : { incoming.first, incoming.second } ) {
    ParticleHistory bhist( beamParticle );
    while( bhist.nextVertex() ) {
      const HepMC::GenVertex* vertex= bhist.getVertex();
      for( HepMC::GenVertex::particles_out_const_iterator
	     poutit= vertex->particles_out_const_begin();
	   poutit != vertex->particles_out_const_end(); poutit++ ) {
	const HepMC::GenParticle* outparticle= *poutit;
	if( outparticle->pdg_id() == 22 and outparticle->status() == 1 ) {
	  ISRphotons.push_back( outparticle );
	}
      }
    }
  }
}

// event_scale is m_inv of sum of beam particle 4-momenta:
Double_t event_scale( const HepMC::GenEvent & event ) {
  std::pair<HepMC::GenParticle*,HepMC::GenParticle*> beamParticles= event.beam_particles();
  const HepMC::FourVector mom1( beamParticles.first->momentum() );
  const HepMC::FourVector mom2( beamParticles.second->momentum() );
  HepMC::FourVector momsum( mom1.px() + mom2.px(),
			    mom1.py() + mom2.py(),
			    mom1.pz() + mom2.pz(),
			    mom1.e() + mom2.e() );
  return momsum.m();
}

// All "undecayed" particles not in ISR photon list are hadron level:
class IsParticleEqual {
  const HepMC::GenParticle* particle;
public:
  IsParticleEqual( const HepMC::GenParticle* p ) : particle(p) {}
  bool operator() ( const HepMC::GenParticle* p ) { return *p == *particle; }
};
void HepMC2Reader::getHadron() {
  vtlv.clear();
  TLorentzVector sum;
  for( HepMC::GenEvent::particle_const_iterator it= event.particles_begin();
       it != event.particles_end(); it++ ) {
    const HepMC::GenParticle* particle= *it;
    if( particle->is_undecayed() and
	std::find_if( ISRphotons.begin(), ISRphotons.end(),
		      IsParticleEqual( particle ) ) == ISRphotons.end() ) {   
      const HepMC::FourVector mom( particle->momentum() );
      TLorentzVector tlv( mom.px(), mom.py(), mom.pz(), mom.e() );
      vtlv.push_back( tlv );
      sum+= tlv;
    }
  }
  for( const HepMC::GenParticle* isrphoton : ISRphotons ) {
    const HepMC::FourVector mom( isrphoton->momentum() );
    TLorentzVector tlv( mom.px(), mom.py(), mom.pz(), mom.e() );
    sum+= tlv;
  }
  if( fabs( sum.Px() ) > 1.0e-6 or fabs( sum.Py() ) > 1.0e-6 or
      fabs( sum.Pz() ) > 1.0e-6 or
      fabs( sum.E() - event_scale( event ) ) > 1.0e-6 ) {
    throw std::runtime_error( "HepMC2Reader::getHadron: 4-vector sum mismatch" );
  }
  return;
}


// Find "hadronisation vertex" and collect partons from "incoming":
void HepMC2Reader::getParton() {
  vtlv.clear();
  TLorentzVector sum;
  for( HepMC::GenEvent::vertex_const_iterator vit= event.vertices_begin();
       vit != event.vertices_end(); vit++ ) {
    const HepMC::GenVertex* vertex= *vit;
    bool allpartonsin= true;
    for( HepMC::GenVertex::particles_in_const_iterator
	   pinit= vertex->particles_in_const_begin();
	 pinit != vertex->particles_in_const_end(); pinit++ ) {
      const HepMC::GenParticle* particle= *pinit;
      allpartonsin= ( allpartonsin and abs( particle->pdg_id() ) < 23 );
    }
    bool allhadronsout= true;
    for( HepMC::GenVertex::particles_out_const_iterator
	   poutit= vertex->particles_out_const_begin();
	 poutit != vertex->particles_out_const_end(); poutit++ ) {
      const HepMC::GenParticle* particle= *poutit;
      allhadronsout= ( allhadronsout and
		       ( abs( particle->pdg_id() ) > 100 or particle->pdg_id() == 22 ) 
		       and ( particle->status() == 1 or particle->status() == 2 ) );
    }
    if( allpartonsin and allhadronsout ) {
      for( HepMC::GenVertex::particles_in_const_iterator
	     pinit= vertex->particles_in_const_begin();
	   pinit != vertex->particles_in_const_end(); pinit++ ) {
	const HepMC::FourVector mom( (*pinit)->momentum() );
	TLorentzVector tlv( mom.px(), mom.py(), mom.pz(), mom.e() );
	vtlv.push_back( tlv );
	sum+= tlv;
      }
      break;
    }
  }
  for( const HepMC::GenParticle* isrphoton : ISRphotons ) {
    const HepMC::FourVector mom( isrphoton->momentum() );
    TLorentzVector tlv( mom.px(), mom.py(), mom.pz(), mom.e() );
    sum+= tlv;
  }
  if( fabs( sum.Px() ) > 1.0e-6 or fabs( sum.Py() ) > 1.0e-6 or
      fabs( sum.Pz() ) > 1.0e-6 or
      fabs( sum.E() - event_scale( event ) ) > 1.0e-6 ) {
    event.print();
    std::cout << "HepMC2Reader::getParton: ISR photons" << std::endl;
    for( const TLorentzVector & isrphoton : vtlv ) isrphoton.Print();
    std::cout << "HepMC2Reader::getParton: partons" << std::endl;
    for( const TLorentzVector & tlv : vtlv ) tlv.Print();
    throw std::runtime_error( "HepMC2Reader::getParton: 4-vector sum mismatch" );
  }
  return;
}

// Get ISR photons:
void HepMC2Reader::getIsr() {
  vtlv.clear();
  for( const HepMC::GenParticle* isrphoton : ISRphotons ) {
    const HepMC::FourVector mom( isrphoton->momentum() );
    TLorentzVector tlv( mom.px(), mom.py(), mom.pz(), mom.e() );
    vtlv.push_back( tlv );
  }
  return;
}

// Return 4-vectors at hadron or parton level:
const std::vector<TLorentzVector>
HepMC2Reader::GetLorentzVectors( const std::string & opt ) {
  if( cacheIsValid.at( opt ) ) return vtlvCache.at( opt );
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

bool HepMC2Reader::GetEvent( int ievent ) {
  return GetNextEvent();
}

// Read next event, return false to stop:
bool HepMC2Reader::GetNextEvent( Int_t maxevt ) {

  // Check maxevt:
  if( maxevt > 0 and maxevt == nevent ) {
    std::cout << "HepMC2Reader::GetNextEvent: event " << nevent
	      << " maximum reached, stop" << std::endl;
    return false;
  }
  
  // Three reads and you're out:
  Int_t itry= 0;
  do {
    try {
      event.read( input );
      if( input.eof() ) {
	std::cout << "HepMC2Reader::GetNextEvent: event " << nevent
		  << " input EOF, stop" << std::endl;
	return false;
      }
      else if( input.bad() ) {
	std::cout << "HepMC2Reader::GetNextEvent: event " << nevent
		  << " input badbit set, skip" << std::endl;
      }
      else if( not event.is_valid() ) {
	std::cout << "HepMC2Reader::GetNextEvent: event " << nevent
		  << " not valid, skip" << std::endl;
      }
      else {
	break;
      }
    }
    catch( const HepMC::IO_Exception & ioe ) {
      std::cout << "HepMC2Reader::GetNextEvent: event " << nevent
		<< " cought HepMC::IO_Exception "
		<< ioe.what() << ", skip" << std::endl;
    }
    itry++;
  }
  while( itry < 3 );
  if( itry == 3 ) {
    std::cout << "HepMC2Reader::GetNextEvent: event " << nevent
	      << " tried reading 3 times, stop" << std::endl;
    return false;
  }

  // Looks good, prepare event processing:
  for( const auto & keyValue : cacheIsValid ) cacheIsValid[keyValue.first]= false;
  findISRphotons();
  if( ( nevent % 1000 ) == 0 and not nevent == 0 ) {
    std::cout << "HepMC2Reader::GetNextEvent: event " << nevent << std::endl;
  }
  nevent++;
  return true;

}

bool HepMC2Reader::MCNonRad() {
  double sumE= 0.0;
  for( const HepMC::GenParticle* p : ISRphotons ) sumE+= p->momentum().e();
  bool result= false;
  if( sumE < 1.0 ) result= true;
  return result;
}

void HepMC2Reader::printParticlesVertices() {
  Int_t ip= 0;
  for( HepMC::GenEvent::particle_const_iterator it= event.particles_begin();
       it != event.particles_end(); it++ ) {
    const HepMC::GenParticle* particle= *it;
    const HepMC::FourVector mom( particle->momentum() );
    std::cout << ip << " " << particle->pdg_id() << " " << particle->status() << " "
	      << particle->barcode()-10001 << " "
	      << particle->generated_mass() << " "
	      << mom.x() << " "
	      << mom.y() << " "
	      << mom.z() << " "
	      << mom.t() << " "
	      << std::endl;
    ip++;
  }
  std::cout << "ISR photons" << std::endl;
  for( const HepMC::GenParticle* photon : ISRphotons ) {
    const HepMC::FourVector mom( photon->momentum() );
    std::cout << photon->barcode()-10001 << " "
	      << photon->pdg_id() << " "
	      << photon->status() << " "
	      << photon->generated_mass() << " "
	      << mom.x() << " "
	      << mom.y() << " "
	      << mom.z() << " "
	      << mom.t()
	      << std::endl;    
  }  
  Int_t iv= 0;
  for( HepMC::GenEvent::vertex_const_iterator vit= event.vertices_begin();
       vit != event.vertices_end(); vit++ ) {
    const HepMC::GenVertex* vertex= *vit;
    std::cout << "vertex " << iv;
    if( vertex->particles_in_size() > 0 ) {
      std::cout << " incoming";
      for( HepMC::GenVertex::particles_in_const_iterator
	     pinit= vertex->particles_in_const_begin();
	   pinit != vertex->particles_in_const_end(); pinit++ ) {
	const HepMC::GenParticle* particle= *pinit;
	std::cout << " " << particle->barcode()-10001;
      }
    }
    if( vertex->particles_out_size() > 0 ) {
      std::cout << " outgoing";
      for( HepMC::GenVertex::particles_out_const_iterator
	     poutit= vertex->particles_out_const_begin();
	   poutit != vertex->particles_out_const_end(); poutit++ ) {
	const HepMC::GenParticle* particle= *poutit;
	std::cout << " " << particle->barcode()-10001;
      }
    }
    std::cout << std::endl;
    iv++;
  }  
}

