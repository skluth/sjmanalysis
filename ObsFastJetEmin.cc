
#include "ObsFastJetEmin.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"

ObsFastJetEmin::ObsFastJetEmin( const string& name, const string& algo, Int_t njet, 
				Double_t rval, 
				const vector<Double_t>& eminfrac, 
				const vector<Analysis>& variations ) :
  ObsJetrate( name ), Algorithm(algo), Jetrate(njet), Rvalue(rval), EminFractions(eminfrac) {
  addAnalyses( variations, EminFractions );
}

void ObsFastJetEmin::fill( NtupleReader* ntr, const Analysis& variation ) {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( variation.getReco() );
  size_t n= EminFractions.size();
  vector<Double_t> NJets( n );
  TFastJet tfj( vtlv, Algorithm.c_str(), Rvalue );
  Double_t Evis= tfj.Evis();
  for( size_t i= 0; i < n; i++ ) {
    vector<TLorentzVector> incljets= tfj.inclusive_jets( EminFractions[i]*Evis );
    NJets[i]= incljets.size();
  }
  getAndFillJetrateDataStructure( NJets, Jetrate, variation.getTag() );
  return;
}

