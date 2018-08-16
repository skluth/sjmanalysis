#ifndef LEPYCUTCALCULATOR_HH
#define LEPYCUTCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>

class LEPYcutCalculator: public JetrateCalculator {
  std::string algorithm;
public:
  LEPYcutCalculator( const std::string& );
  ~LEPYcutCalculator() {}
  std::vector<Double_t> getValues( NtupleReader*, 
				   const std::vector<Double_t> &,
				   const std::string & ) const;
  void print() const;
};

#endif
