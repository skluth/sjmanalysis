
#include "YnmCalculator.hh"
#include <string>
#include "NtupleReader.hh"

#include <iostream>

Double_t YnmCalculator::getValue( NtupleReader* ntr, 
				  const std::string& reco ) const {
  return ntr->getYmerge( algorithm, reco, njet );
}
