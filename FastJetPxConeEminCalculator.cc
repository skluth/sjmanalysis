
#include "FastJetPxConeEminCalculator.hh"

#include "NtupleReader.hh"
#include "TFastJet.hh"

#include "TLorentzVector.h"

#include <iostream>

using std::string;
using std::vector;

FastJetPxConeEminCalculator::FastJetPxConeEminCalculator( const Double_t R ) :
  RValue(R) {}

void FastJetPxConeEminCalculator::print() const {
  std::cout << "FastJetPxConeEminCalculator R value: " << RValue << std::endl;
}

vector<Double_t> 
FastJetPxConeEminCalculator::getValues( NtupleReader* ntr, 
					const vector<Double_t>& EminPoints,
					const string& reco ) const {
  const vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( reco );
  size_t n= EminPoints.size();
  vector<Double_t> NJets( n );
  for( size_t i= 0; i < n; i++ ) {
    TFastJet tfj( vtlv, "pxcone", RValue, EminPoints[i] );
    vector<TLorentzVector> incljets= tfj.inclusive_eejets( EminPoints[i] );      
    NJets[i]= incljets.size();
  }
  return NJets;
}

