#ifndef OBSFASTJETEMIN_HH
#define OBSFASTJETEMIN_HH

#include "ObsJetrate.hh"
#include "Analysis.hh"

#include <vector>
using std::vector;
#include <string>
using std::string;

class NtupleReader;

class ObsFastJetEmin : public ObsJetrate {

public:

  ObsFastJetEmin( const string& name, const string& algo, Double_t rval, 
		  const vector<Double_t>& eminfractions, 
		  const vector<Analysis>& variations );
  ~ObsFastJetEmin() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation );
  
private:
  
  string Algorithm;
  Double_t Rvalue;

};


#endif
