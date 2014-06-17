#ifndef OBSDURHAMYMERGE23_HH
#define OBSDURHAMYMERGE23_HH

#include "Rtypes.h"
#include "ObsDifferential.hh"
#include "Analysis.hh"

#include <vector>
#include <string>
using std::vector;
using std::string;

class NtupleReader;

class ObsDurhamYmerge23 : public ObsDifferential {

public:

  ObsDurhamYmerge23( const vector<Double_t>& bins, 
		     const vector<Analysis>& variations );
  ~ObsDurhamYmerge23() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation );

};

#endif
