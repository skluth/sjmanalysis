#ifndef PXTHRUSTCALCULATOR_HH
#define PXTHRUSTCALCULATOR_HH

#include "DifferentialCalculator.hh"

class NtupleReader;
#include <string>

class PxThrustCalculator: public DifferentialCalculator {
public:
  PxThrustCalculator() {}
  ~PxThrustCalculator() {}
  Double_t getValue( NtupleReader*, const std::string& ) const;
};

#endif
