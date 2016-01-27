
#include "FastJetYcutCalculator.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
#include "TMath.h"
#include "TLorentzVector.h"

FastJetYcutCalculator::FastJetYcutCalculator( const string& algo ) : 
  algorithm(algo) {}

vector<Double_t> 
FastJetYcutCalculator::getValues( NtupleReader* ntr, 
				  const vector<Double_t>& Ycutpoints,
				  const string& reco) const {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( reco );
  size_t n= Ycutpoints.size();
  vector<Double_t> NJets( n );
  TFastJet tfj( vtlv, algorithm.c_str() );
  for( size_t i= 0; i < n; i++ ) {
    NJets[i]= tfj.njets( TMath::Power( 10.0, -Ycutpoints[i] ) );
  }
  return NJets;
}

