#ifndef FASTJETPXCONEEMINCALCULATOR_HH
#define FASTJETPXCONEEMINCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>
using std::string;
using std::vector;

class FastJetPxConeEminCalculator: public JetrateCalculator {
  Double_t RValue;
public:
  FastJetPxConeEminCalculator( Double_t R );
  ~FastJetPxConeEminCalculator() {}
  vector<Double_t> getValues( NtupleReader*, 
			      const vector<Double_t>&,
			      const string& ) const;
};

#endif
