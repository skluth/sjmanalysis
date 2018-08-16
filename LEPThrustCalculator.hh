#ifndef LEPTHRUSTCALCULATOR_HH
#define LEPTHRUSTCALCULATOR_HH

#include "DifferentialCalculator.hh"

class NtupleReader;
#include <string>

class LEPThrustCalculator: public DifferentialCalculator {
public:
  LEPThrustCalculator() {}
  ~LEPThrustCalculator() {}
  Double_t getValue( NtupleReader*, const std::string& ) const;
};

#endif
