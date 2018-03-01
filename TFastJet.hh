#ifndef TFASTJET_HH
#define TFASTJET_HH

#include "TLorentzVector.h"
#include <vector>
using std::vector;

namespace fastjet {
  class ClusterSequence;
  class PseudoJet;
}

class TFastJet {

public:

  TFastJet() {}

  TFastJet( const vector<TLorentzVector>&, const char* jetalg="antikt",
	    const double& R=0.4, const vector<int>* vindx= 0,
	    const double Emin=0.0 );

  virtual ~TFastJet();

  const vector<TLorentzVector> inclusive_jets( const double& ptmin=0.0 );
  const vector<TLorentzVector> inclusive_eejets( const double& Emin=0.0 );
  const vector<TLorentzVector> exclusive_jets( const int njets );
  double ymerge( int );
  int njets( double );
  double Evis();

private:

  const vector<TLorentzVector>
  copyPseudoJetsToLorentzVectors( const vector<fastjet::PseudoJet>&,
				  const double Emin=0.0 );

  vector< vector<int> > constituents( const vector<fastjet::PseudoJet>& pjets );

  fastjet::ClusterSequence* clusseq;

};


#endif
