
#include "LEPThrustCalculator.hh"
#include <string>
#include <stdexcept>
#include "NtupleReader.hh"
#include "LEPNtupleReader.hh"

Double_t LEPThrustCalculator::getValue( NtupleReader* ntr, 
					const std::string& reco ) const {
  LEPNtupleReader* lepntr= dynamic_cast<LEPNtupleReader*>( ntr );
  if( lepntr == nullptr ) {
    throw std::runtime_error( "LEPThrustCalculator::getValue: no LEP ntuple" );
  }
  return lepntr->getThrust( reco );
}
