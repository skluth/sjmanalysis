#ifndef FASTJETEMINCALCULATOR_HH
#define FASTJETEMINCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>

class FastJetEminCalculator: public JetrateCalculator {
  std::string algorithm;
  Double_t Rvalue;
public:
  FastJetEminCalculator( const std::string &, Double_t );
  ~FastJetEminCalculator() {}
  std::vector<Double_t> getValues( NtupleReader*, 
				   const std::vector<Double_t> &,
				   const std::string & ) const;
  void print() const;
};

#endif
