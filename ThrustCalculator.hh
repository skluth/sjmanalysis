#ifndef THRUSTCALCULATOR_HH
#define THRUSTCALCULATOR_HH

#include "DifferentialCalculator.hh"

class NtupleReader;
#include <string>

class ThrustCalculator: public DifferentialCalculator {
public:
  ThrustCalculator() {}
  ~ThrustCalculator() {}
  Double_t getValue( NtupleReader*, const std::string& ) const;
};

#endif
