#ifndef OBSJETRATE_HH
#define OBSJETRATE_HH

#include "Observable.hh"
#include "Rtypes.h"
#include "Analysis.hh"
#include <vector>
using std::vector;
#include <string>
using std::string;

class NtupleReader;

class ObsJetrate : public Observable {

public:

  ObsJetrate( string name );
  ObsJetrate() {}
  ~ObsJetrate() {}
  void addAnalyses( const vector<Analysis>& variations, const vector<Double_t>& points );
  virtual void fill( NtupleReader* ntr, const Analysis& variation ) = 0;

protected:

  void getAndFillJetrateDataStructure( vector<Double_t> NJets, Int_t Jetrate,
				       const string& tag );

};

#endif
