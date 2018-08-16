
#include "LEPYnmCalculator.hh"
#include <string>
#include "NtupleReader.hh"
#include "LEPNtupleReader.hh"

#include <iostream>
#include <stdexcept>

Double_t LEPYnmCalculator::getValue( NtupleReader* ntr, 
				  const std::string& reco ) const {
  LEPNtupleReader* lepntr= dynamic_cast<LEPNtupleReader*>( ntr );
  if( lepntr == nullptr ) {
    throw std::runtime_error( "LEPYnmCalculator::getValue: no LEP ntuple" );
  }
  return lepntr->getYmerge( algorithm, reco, njet );
}
