#ifndef JETRATECALCULATOR_HH
#define JETRATECALCULATOR_HH

#include "Rtypes.h"
#include <string>
#include <vector>
using std::string;
using std::vector;

class NtupleReader;

class JetrateCalculator {
public:
  vector<Double_t> virtual getValues( NtupleReader*, 
				      const vector<Double_t>&,
				      const string& ) const = 0;
};

#endif
