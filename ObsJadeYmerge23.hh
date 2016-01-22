#ifndef OBSJADEYMERGE23_HH
#define OBSJADEYMERGE23_HH

#include "Rtypes.h"
#include "ObsDifferential.hh"
#include "Analysis.hh"

#include <vector>
using std::vector;

class NtupleReader;

class ObsJadeYmerge23 : public ObsDifferential {

public:

  ObsJadeYmerge23( const vector<Double_t>& bins, 
		   const vector<Analysis>& variations );
  ~ObsJadeYmerge23() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation );

};

#endif
