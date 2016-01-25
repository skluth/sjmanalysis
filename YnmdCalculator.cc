
#include "YnmdCalculator.hh"
#include <string>
#include "NtupleReader.hh"

#include <iostream>

Double_t YnmdCalculator::getValue( NtupleReader* ntr, 
				   const std::string& reco ) const {
  return -TMath::Log10( ntr->getYmergeD( reco, njet ) );
}
