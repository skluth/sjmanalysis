
#include "FastJetPxConeRCalculator.hh"

#include "NtupleReader.hh"
#include "TFastJet.hh"

#include "TLorentzVector.h"

#include <iostream>

using std::string;
using std::vector;

FastJetPxConeRCalculator::FastJetPxConeRCalculator( const Double_t Emin ) :
  EminValue(Emin) {}

void FastJetPxConeRCalculator::print() const {
  std::cout << "FastJetPxConeRCalculator Emin value: "
	    << EminValue << " GeV" << std::endl;
}

vector<Double_t> 
FastJetPxConeRCalculator::getValues( NtupleReader* ntr, 
				     const vector<Double_t>& RPoints,
				     const string& reco ) const {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( reco );
  size_t n= RPoints.size();
  vector<Double_t> NJets( n );
  for( size_t i= 0; i < n; i++ ) {
    TFastJet tfj( vtlv, "pxcone", RPoints[i], EminValue );
    vector<TLorentzVector> incljets= tfj.inclusive_eejets( EminValue );
    NJets[i]= incljets.size();
  }
  return NJets;
}

