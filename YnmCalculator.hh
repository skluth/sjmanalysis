#ifndef YNMCALCULATOR_HH
#define YNMCALCULATOR_HH

#include "DifferentialCalculator.hh"

class NtupleReader;
#include <string>

class YnmCalculator: public DifferentialCalculator {
  std::string algorithm;
  Int_t njet;
public:
  YnmCalculator( const std::string & alg, Int_t n ) : algorithm(alg), njet(n) {}
  ~YnmCalculator() {}
  Double_t getValue( NtupleReader*, const std::string& ) const;
};

#endif
