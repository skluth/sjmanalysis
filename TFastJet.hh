#ifndef TFASTJET_HH
#define TFASTJET_HH

#include <vector>
#include "TLorentzVector.h"
using std::vector;

class TParticle;

namespace fastjet {
  class ClusterSequence;
  class PseudoJet;
}

// rootcint can't cope with fastjet's nested base classes, so have to
// fool it here:
#ifdef __CINT__
namespace fastjet {
  namespace JetDefinition {
    class Plugin;
  }
}
#else
#include "fastjet/JetDefinition.hh"
#endif

class TFastJet {

public:

  TFastJet() {}

  TFastJet( const vector<TParticle>& );
  TFastJet( const vector<TLorentzVector>&, const char* jetalg="antikt",
	    const double& R=0.4, const vector<int>* vindx= 0 );

  virtual ~TFastJet();

  vector<TLorentzVector>& inclusive_jets( const double& ptmin=0.0 );
  vector<TLorentzVector>& exclusive_jets( const int njets );
  vector< vector<int> >& constituents();
  double ymerge( int );
  int njets( double );
  double Evis();

private:

  vector<TLorentzVector>& copyPseudoJetsToLorentzVectors();
  fastjet::ClusterSequence* clusseq;
  fastjet::JetDefinition::Plugin* plugin;
  vector<fastjet::PseudoJet>* pjets;

};

#endif
