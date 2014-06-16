#ifndef OBSFASTJETR_HH
#define OBSFASTJETR_HH

#include "ObsJetrate.hh"
#include "Analysis.hh"

#include <vector>
using std::vector;
#include <string>
using std::string;

class NtupleReader;

class ObsFastJetR : public ObsJetrate {

public:

  ObsFastJetR( const string& name, const string& algo, Int_t njet, Double_t eminfrac, 
	       const vector<Double_t>& rvals, 
	       const vector<Analysis>& variations );
  ~ObsFastJetR() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation );

private:

  string Algorithm;
  Int_t Jetrate;
  Double_t EminFraction;
  vector<Double_t> Rvalues;

};


#endif
