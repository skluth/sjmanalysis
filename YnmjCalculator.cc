
#include "YnmjCalculator.hh"
#include <string>
#include "NtupleReader.hh"

#include <iostream>

Double_t YnmjCalculator::getValue( NtupleReader* ntr, 
				   const std::string& reco ) const {
  return -TMath::Log10( ntr->getYmergeE( reco, njet ) );
}
