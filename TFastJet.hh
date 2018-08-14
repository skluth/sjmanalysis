#ifndef TFASTJET_HH
#define TFASTJET_HH

#include "TLorentzVector.h"
#include <vector>

namespace fastjet {
  class ClusterSequence;
  class PseudoJet;
}

class JetDefinition_Recombiner;
class JetDefinition_Plugin;

class TFastJet {

public:

  TFastJet() {}

  TFastJet( const std::vector<TLorentzVector>&,
	    const char* jetalg="antikt",
	    const double R=0.4,
	    const double Emin=0.0 );

  virtual ~TFastJet();

  const std::vector<TLorentzVector> inclusive_jets( const double ptmin=0.0 );
  const std::vector<TLorentzVector> inclusive_eejets( const double Emin=0.0 );
  const std::vector<TLorentzVector> exclusive_eejets( const int njets );
  const std::vector<TLorentzVector>
  inclusive_jets( std::vector< std::vector<int> > & vindx,
		  const double ptmin=0.0 );
  const std::vector<TLorentzVector>
  inclusive_eejets( std::vector< std::vector<int> > & vindx,
		    const double Emin=0.0 );
  const std::vector<TLorentzVector>
  exclusive_eejets( std::vector< std::vector<int> > & vindx,
		    const int njets );
  double ymerge( int );
  int njets( double );
  double Evis();

private:

  const std::vector<TLorentzVector>
  copyPseudoJetsToLorentzVectors( const std::vector<fastjet::PseudoJet>&,
				  const double Emin=0.0 );

  std::vector< std::vector<int> >
  constituents( const std::vector<fastjet::PseudoJet>& pjet );

  fastjet::ClusterSequence* clusseq;

};

#endif
