#ifndef FASTJETYCUTCALCULATOR_HH
#define FASTJETYCUTCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>
using std::string;
using std::vector;

class FastJetYcutCalculator: public JetrateCalculator {
  string algorithm;
public:
  FastJetYcutCalculator( const string& );
  ~FastJetYcutCalculator() {}
  vector<Double_t> getValues( NtupleReader*, 
			      const vector<Double_t>&,
			      const string& ) const;
  void print() const;
};

#endif
