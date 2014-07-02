#ifndef OBSMR_HH
#define OBSMR_HH

#include "Rtypes.h"
#include "ObsDifferential.hh"
#include "Analysis.hh"

#include <vector>
using std::vector;

class NtupleReader;

class ObsMr : public ObsDifferential {

public:

  ObsMr( const vector<Double_t>& bins, 
	 const vector<Analysis>& variations,
	 Double_t y34c=0.0045, Double_t y34y23c=0.5 );
  ~ObsMr() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation );

private:

  Double_t getMrvalue( NtupleReader* ntr, const string& reco );

  Double_t y34cut;
  Double_t y34y23cut;

};

#endif
