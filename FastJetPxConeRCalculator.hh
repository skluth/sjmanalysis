#ifndef FASTJETPXCONERCALCULATOR_HH
#define FASTJETPXCONERCALCULATOR_HH

#include "JetrateCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>

class FastJetPxConeRCalculator: public JetrateCalculator {
  Double_t EminValue;
public:
  FastJetPxConeRCalculator( const Double_t Emin );
  ~FastJetPxConeRCalculator() {}
  std::vector<Double_t> getValues( NtupleReader*, 
				   const std::vector<Double_t> &,
				   const std::string & ) const;
  void print() const;
};

#endif
