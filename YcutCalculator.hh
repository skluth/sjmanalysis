#ifndef YCUTCALCULATOR_HH
#define YCUTCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>

class YcutCalculator: public JetrateCalculator {
  std::string algorithm;
public:
  YcutCalculator( const std::string& );
  ~YcutCalculator() {}
  std::vector<Double_t> getValues( NtupleReader*, 
				   const std::vector<Double_t> &,
				   const std::string & ) const;
  void print() const;
};

#endif
