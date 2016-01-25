#ifndef YNMDCALCULATOR_HH
#define YNMDCALCULATOR_HH

#include "DifferentialCalculator.hh"

class NtupleReader;
#include <string>

class YnmdCalculator: public DifferentialCalculator {
  Int_t njet;
public:
  YnmdCalculator( Int_t n ) : njet(n) {}
  ~YnmdCalculator() {}
  Double_t getValue( NtupleReader*, const std::string& ) const;
};

#endif
