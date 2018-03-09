#ifndef FASTJETPXCONEEMINCALCULATOR_HH
#define FASTJETPXCONEEMINCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>

class FastJetPxConeEminCalculator: public JetrateCalculator {
  Double_t RValue;
public:
  FastJetPxConeEminCalculator( const Double_t R );
  ~FastJetPxConeEminCalculator() {}
  std::vector<Double_t> getValues( NtupleReader*, 
				   const std::vector<Double_t> &,
				   const std::string & ) const;
  void print() const;
};

#endif
