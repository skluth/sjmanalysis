#ifndef FASTJETYCUTCALCULATOR_HH
#define FASTJETYCUTCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>

class FastJetYcutCalculator: public JetrateCalculator {
  std::string algorithm;
public:
  FastJetYcutCalculator( const std::string & );
  ~FastJetYcutCalculator() {}
  std::vector<Double_t> getValues( NtupleReader*, 
				   const std::vector<Double_t> &,
				   const std::string & ) const;
  void print() const;
};

#endif
