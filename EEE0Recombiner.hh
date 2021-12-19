#ifndef EEE0RECOMBINER_HH
#define EEE0RECOMBINER_HH

#include "fastjet/JetDefinition.hh"
#include <string>

class EEE0Recombiner: public fastjet::JetDefinition::Recombiner {
public:
  EEE0Recombiner() {
    std::cout << "EEE0Recombiner::EEE0Recombiner: created" << std::endl;
  }
  ~EEE0Recombiner() {}
  virtual std::string description() const{
    return "E0 scheme for EE";
  }
  virtual void recombine( const fastjet::PseudoJet& pa,
                          const fastjet::PseudoJet& pb, 
                          fastjet::PseudoJet& pab ) const {
    pab.reset( pa.px() + pb.px(), pa.py() + pb.py(),
               pa.pz() + pb.pz(), pa.E() + pb.E() );
    // double rescale= pab.E()/TMath::Sqrt( pab.perp2() + pab.pz()*pab.pz() );
    // pab.reset_momentum( rescale*pab.px(), rescale*pab.py(), rescale*pab.pz(), pab.E() );
    return;
  }
};

#endif
