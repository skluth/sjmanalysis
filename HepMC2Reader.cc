
#include "HepMC2Reader.hh"

#include "HepMC/GenParticle.h"
#include "HepMC/SimpleVector.h"

HepMC2Reader::HepMC2Reader( const std::string& filename ) :
  input( filename ) {}

bool HepMC2Reader::GetEvent( Int_t ievnt ) {
  event.read( input );
  return event.is_valid();
}

const std::vector<TLorentzVector> HepMC2Reader::GetLorentzVectors( const std::string & opt ) {
  std::vector<TLorentzVector> vtlv;
  for( HepMC::GenEvent::particle_const_iterator it= event.particles_begin();
       it != event.particles_end(); it++ ) {
    const HepMC::GenParticle* particle= *it;
    if( particle->is_undecayed() ) {
      const HepMC::FourVector mom( particle->momentum() );
      TLorentzVector tlv( mom.px(), mom.py(), mom.pz(), mom.e() );
      vtlv.push_back( tlv );
    }
  }
  return vtlv;
}
