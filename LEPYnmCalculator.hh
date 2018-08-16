#ifndef LEPYNMCALCULATOR_HH
#define LEPYNMCALCULATOR_HH

#include "DifferentialCalculator.hh"

class NtupleReader;
#include <string>

class LEPYnmCalculator: public DifferentialCalculator {
  std::string algorithm;
  Int_t njet;
public:
  LEPYnmCalculator( const std::string & alg, Int_t n ) : algorithm(alg), njet(n) {}
  ~LEPYnmCalculator() {}
  Double_t getValue( NtupleReader*, const std::string& ) const;
};

#endif
