#ifndef OBSTHRUST_HH
#define OBSTHRUST_HH

#include "Rtypes.h"
#include "ObsDifferential.hh"
#include "Analysis.hh"

#include <vector>
#include <string>
using std::vector;
using std::string;

class NtupleReader;

class ObsThrust : public ObsDifferential {

public:

  ObsThrust( const vector<Double_t>& bins, 
	     const vector<Analysis>& variations );
  ~ObsThrust() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation );

};

#endif
