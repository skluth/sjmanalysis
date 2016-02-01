
#include "FastJetEminCalculator.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
#include "TLorentzVector.h"

FastJetEminCalculator::FastJetEminCalculator( const string& algo, 
					      Double_t R) : 
  algorithm(algo), Rvalue(R) {}

vector<Double_t> 
FastJetEminCalculator::getValues( NtupleReader* ntr, 
				  const vector<Double_t>& Eminfpoints,
				  const string& reco ) const {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( reco );
  size_t n= Eminfpoints.size();
  vector<Double_t> NJets( n );
  TFastJet tfj( vtlv, algorithm.c_str(), Rvalue );
  Double_t Evis= tfj.Evis();
  for( size_t i= 0; i < n; i++ ) {
    // vector<TLorentzVector> incljets= tfj.inclusive_jets( Eminfpoints[i]*Evis );
    vector<TLorentzVector> incljets= tfj.inclusive_eejets( Eminfpoints[i]*Evis );
    NJets[i]= incljets.size();
  }
  return NJets;
}

