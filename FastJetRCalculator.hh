#ifndef FASTJETRCALCULATOR_HH
#define FASTJETRCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>
using std::string;
using std::vector;

class FastJetRCalculator: public JetrateCalculator {
  //vector<Double_t> points;
  string algorithm;
  Double_t EminFraction;
public:
  FastJetRCalculator( const string&, Double_t );
  ~FastJetRCalculator() {}
  vector<Double_t> getValues( NtupleReader*, 
			      const vector<Double_t>&,
			      const string& ) const;
};

#endif
