#ifndef TFASTJET_HH
#define TFASTJET_HH

#include "TLorentzVector.h"
#include <vector>
using std::vector;
#include <string>
using std::string;

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
    class Recombiner;
  }
}
#else
#include "fastjet/JetDefinition.hh"
#endif

class TFastJet {

public:

  TFastJet() {}

  //  TFastJet( const vector<TParticle>& );
  TFastJet( const vector<TLorentzVector>&, const char* jetalg="antikt",
	    const double& R=0.4, const vector<int>* vindx= 0 );

  virtual ~TFastJet();

  const vector<TLorentzVector>& inclusive_jets( const double& ptmin=0.0 );
  const vector<TLorentzVector>& exclusive_jets( const int njets );
  vector< vector<int> >& constituents();
  double ymerge( int );
  int njets( double );
  double Evis();

private:

  const vector<TLorentzVector>& copyPseudoJetsToLorentzVectors();

  fastjet::ClusterSequence* clusseq;
  fastjet::JetDefinition::Plugin* plugin;
  vector<fastjet::PseudoJet>* pjets;

};

class EEE0Recombiner: public fastjet::JetDefinition::Recombiner {

public:

  EEE0Recombiner() {}
  ~EEE0Recombiner() {}
  virtual string description() const;
  virtual void recombine( const fastjet::PseudoJet& pa,
			  const fastjet::PseudoJet& pb, 
			  fastjet::PseudoJet& pab ) const;

};

#endif
