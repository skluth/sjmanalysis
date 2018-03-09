#ifndef FASTJETRCALCULATOR_HH
#define FASTJETRCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>

class FastJetRCalculator: public JetrateCalculator {
  std::string algorithm;
  Double_t EminFraction;
public:
  FastJetRCalculator( const std::string &, Double_t );
  ~FastJetRCalculator() {}
  std::vector<Double_t> getValues( NtupleReader*, 
				   const std::vector<Double_t> &,
				   const std::string & ) const;
  void print() const;
};

#endif
