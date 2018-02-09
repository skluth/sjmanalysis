#ifndef YCUTCALCULATOR_HH
#define YCUTCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>
using std::string;
using std::vector;

class YcutCalculator: public JetrateCalculator {
  string algorithm;
public:
  YcutCalculator( const string& );
  ~YcutCalculator() {}
  vector<Double_t> getValues( NtupleReader*, 
			      const vector<Double_t>&,
			      const string& ) const;
  void print() const;
};

#endif
