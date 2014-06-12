#ifndef OBSFASTJETDIFF_HH
#define OBSFASTJETDIFF_HH

#include "Rtypes.h"
#include "ObsDifferential.hh"
#include "Analysis.hh"

#include <vector>
#include <string>
using std::vector;
using std::string;
using std::map;

class NtupleReader;

class ObsFastJetDiff : public ObsDifferential {

public:

  ObsFastJetDiff( const string& name, const string& algo, Int_t njet,  
		  const vector<Double_t>& bins, 
		  const vector<Analysis>& variations );
  ~ObsFastJetDiff() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation );

private:
  
  string Algorithm;
  Int_t Jetlower;

};

#endif
