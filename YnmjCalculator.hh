#ifndef YNMjCALCULATOR_HH
#define YNMjCALCULATOR_HH

#include "DifferentialCalculator.hh"

class NtupleReader;
#include <string>

class YnmjCalculator: public DifferentialCalculator {
  Int_t njet;
public:
  YnmjCalculator( Int_t n ) : njet(n) {}
  ~YnmjCalculator() {}
  Double_t getValue( NtupleReader*, const std::string& ) const;
};

#endif
