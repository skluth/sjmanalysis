
#include "ObsFastJetYcut.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
#include "TMath.h"

ObsFastJetYcut::ObsFastJetYcut( const string& name, const string& algo, Int_t njet, 
				const vector<Double_t>& ycuts, 
				const vector<Analysis>& variations ) :
  ObsJetrate( name ), Algorithm(algo), Jetrate(njet), YcutPoints(ycuts) {
  addAnalyses( variations, YcutPoints );
}

void ObsFastJetYcut::fill( NtupleReader* ntr, const Analysis& variation ) {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( variation.getReco() );
  size_t n= YcutPoints.size();
  vector<Double_t> NJets( n );
  TFastJet tfj( vtlv, Algorithm.c_str() );
  for( size_t i= 0; i < n; i++ ) {
    NJets[i]= tfj.njets( TMath::Power( 10.0, -YcutPoints[i] ) );
  }
  getAndFillJetrateDataStructure( NJets, Jetrate, variation.getTag() );
  return;
}

