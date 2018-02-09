#ifndef FASTJETEMINCALCULATOR_HH
#define FASTJETEMINCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>
using std::string;
using std::vector;

class FastJetEminCalculator: public JetrateCalculator {
  string algorithm;
  Double_t Rvalue;
public:
  FastJetEminCalculator( const string&, Double_t );
  ~FastJetEminCalculator() {}
  vector<Double_t> getValues( NtupleReader*, 
			      const vector<Double_t>&,
			      const string& ) const;
  void print() const;
};

#endif
