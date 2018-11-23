
#include "HepMC2Reader.hh"

#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"
#include <iostream>
#include <vector>
#include <algorithm>


HepMC2Reader::HepMC2Reader( const std::string & filename ) :
  input( filename ) {
  std::cout << "HepMC2Reader::HepMC2Reader: opened hepmc2 file " << filename << std::endl;
}

// Read next event:
bool HepMC2Reader::GetEvent( Int_t ievnt ) {
  event.read( input );
  if( event.is_valid() ) {
    for( const auto & keyValue : cacheIsValid ) cacheIsValid[keyValue.first]= false;
    //  event.print();
    findISRphotons();
    return true;
  }
  return false;
}

bool HepMC2Reader::MCNonRad() {
  double sumE= 0.0;
  for( const HepMC::GenParticle* p : ISRphotons ) sumE+= p->momentum().e();
  bool result= false;
  if( sumE < 1.0 ) result= true;
  return result;
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
  const HepMC::GenParticle* getParticle() { return particle; }
  const HepMC::GenVertex* getVertex() { return particle->end_vertex(); }
};
// Find ISR from beams:
void HepMC2Reader::findISRphotons() {
  std::pair<HepMC::GenParticle*,HepMC::GenParticle*> incoming= event.beam_particles();
  std::vector<ParticleHistory> bhists;
  bhists.push_back( ParticleHistory( incoming.first ) );
  bhists.push_back( ParticleHistory( incoming.second ) );
  ISRphotons.clear();
  for( ParticleHistory & bhist : bhists ) {
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

// Predicate function object for particle matching with find_if
class IsParticleEqual {
  const HepMC::GenParticle* particle;
public:
  IsParticleEqual( const HepMC::GenParticle* p ) : particle(p) {}
  bool operator() ( const HepMC::GenParticle* p ) { return *p == *particle; }
};

// All "undecayed" particles not in ISR photon list are hadron level:
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
  if( fabs( sum.Px() ) > 1.0e-6 or fabs( sum.Py() ) > 1.0e-6 or fabs( sum.Pz() ) > 1.0e-6 or
      fabs( sum.E() - event.event_scale() ) > 1.0e06 ) {
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
  if( fabs( sum.Px() ) > 1.0e-6 or fabs( sum.Py() ) > 1.0e-6 or fabs( sum.Pz() ) > 1.0e-6 or
      fabs( sum.E() - event.event_scale() ) > 1.0e06 ) {

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
const std::vector<TLorentzVector> HepMC2Reader::GetLorentzVectors( const std::string & opt ) {
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

