
#include "FastJetEminCalculator.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
#include "TLorentzVector.h"
#include <iostream>

FastJetEminCalculator::FastJetEminCalculator( const string& algo, 
					      Double_t R ) : 
  algorithm(algo), Rvalue(R) {}

void FastJetEminCalculator::print() const {
  std::cout << "FastJetEminCalculator algorithm: " << algorithm 
	    << ", R value: " << Rvalue << std::endl;
}

vector<Double_t> 
FastJetEminCalculator::getValues( NtupleReader* ntr, 
				  const vector<Double_t>& Eminfpoints,
				  const string& reco ) const {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( reco );
  Double_t Evis= ntr->Evis( vtlv );
  size_t n= Eminfpoints.size();
  vector<Double_t> NJets( n );
  if( algorithm == "pxcone" ) {
    for( size_t i= 0; i < n; i++ ) {
      TFastJet tfj( vtlv, algorithm.c_str(), Rvalue, Eminfpoints[i]*Evis );
      vector<TLorentzVector> incljets= tfj.inclusive_eejets( Eminfpoints[i]*Evis );
      NJets[i]= incljets.size();
    }
  }
  else {
    TFastJet tfj( vtlv, algorithm.c_str(), Rvalue );
    for( size_t i= 0; i < n; i++ ) {
      vector<TLorentzVector> incljets= tfj.inclusive_eejets( Eminfpoints[i]*Evis );
      NJets[i]= incljets.size();
    }
  }
  return NJets;
}

