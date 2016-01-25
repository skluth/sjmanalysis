
#include "ThrustCalculator.hh"
#include <string>
#include "NtupleReader.hh"

Double_t ThrustCalculator::getValue( NtupleReader* ntr, 
				     const std::string& reco ) const {
  return ntr->getThrust( reco );
}
