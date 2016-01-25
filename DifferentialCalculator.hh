#ifndef DIFFERENTIALCALCULATOR_HH
#define DIFFERENTIALCALCULATOR_HH

class NtupleReader;
#include<string>
#include "Rtypes.h"

class DifferentialCalculator {
public:
  Double_t virtual getValue( NtupleReader*, const std::string& ) const = 0;
};

#endif
