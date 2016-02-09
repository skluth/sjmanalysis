#ifndef FASTJETPXCONERCALCULATOR_HH
#define FASTJETPXCONERCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>
using std::string;
using std::vector;

class FastJetPxConeRCalculator: public JetrateCalculator {
  Double_t EminValue;
public:
  FastJetPxConeRCalculator( const Double_t Emin );
  ~FastJetPxConeRCalculator() {}
  vector<Double_t> getValues( NtupleReader*, 
			      const vector<Double_t>&,
			      const string& ) const;
};

#endif
