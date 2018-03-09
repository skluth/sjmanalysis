#ifndef JETRATECALCULATOR_HH
#define JETRATECALCULATOR_HH

#include "Rtypes.h"

#include <string>
#include <vector>

class NtupleReader;

class JetrateCalculator {
public:
  std::vector<Double_t> virtual getValues( NtupleReader*, 
					   const std::vector<Double_t> &,
					   const std::string & ) const = 0;
  virtual void print() const = 0;
};

#endif
