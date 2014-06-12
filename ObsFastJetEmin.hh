#ifndef OBSFASTJETEMIN_HH
#define OBSFASTJETEMIN_HH

#include "ObsJetrate.hh"
#include "Analysis.hh"

#include <vector>
#include <string>
using std::vector;
using std::string;

class NtupleReader;

class ObsFastJetEmin : public ObsJetrate {

public:

  ObsFastJetEmin( const string& name, const string& algo, Int_t njet, Double_t rval, 
		  const vector<Double_t>& eminfrac, 
		  const vector<Analysis>& variations );
  ~ObsFastJetEmin() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation );
  
private:
  
  string Algorithm;
  Int_t Jetrate;
  Double_t Rvalue;
  vector<Double_t> EminFractions;

};


#endif
