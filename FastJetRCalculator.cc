
#include "FastJetRCalculator.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"

#include "TLorentzVector.h"

#include <iostream>

FastJetRCalculator::FastJetRCalculator( const string& algo, 
					Double_t eminf ) : 
  algorithm(algo), EminFraction(eminf) {}

void FastJetRCalculator::print() const {
  std::cout << "FastJetRCalculator algorithm: " << algorithm 
	    << ", Emin/Evis value: " << EminFraction << std::endl;
}

vector<Double_t> 
FastJetRCalculator::getValues( NtupleReader* ntr, 
			       const vector<Double_t>& Rpoints, 
			       const string& reco ) const {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( reco );

  Double_t Evis= ntr->Evis( vtlv );

  size_t n= Rpoints.size();
  vector<Double_t> NJets( n );
  for( size_t i= 0; i < n; i++ ) {
    // TFastJet tfj( vtlv, algorithm.c_str(), Rpoints[i] );
    TFastJet tfj( vtlv, algorithm.c_str(), Rpoints[i], 0, EminFraction*Evis );
    // Double_t Evis= tfj.Evis();
    vector<TLorentzVector> incljets= tfj.inclusive_eejets( EminFraction*Evis );
    NJets[i]= incljets.size();
  }
  return NJets;
}

