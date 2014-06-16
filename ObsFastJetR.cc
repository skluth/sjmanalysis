
#include "ObsFastJetR.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"

ObsFastJetR::ObsFastJetR( const string& name, const string& algo, Int_t njet, Double_t eminfrac,
			  const vector<Double_t>& rvals, 
			  const vector<Analysis>& variations ) :
  ObsJetrate( name ), Algorithm(algo), Jetrate(njet), EminFraction(eminfrac), 
  Rvalues(rvals) {
  addAnalyses( variations, Rvalues );
}

void ObsFastJetR::fill( NtupleReader* ntr, const Analysis& variation ) {
  vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( variation.getReco() );
  size_t n= Rvalues.size();
  vector<Double_t> NJets( n );
  for( size_t i= 0; i < n; i++ ) {
    TFastJet tfjakt( vtlv, Algorithm.c_str(), Rvalues[i] );
    Double_t Evis= tfjakt.Evis();
    vector<TLorentzVector> incljets= tfjakt.inclusive_jets( EminFraction*Evis );
    NJets[i]= incljets.size();
  }
  getAndFillJetrateDataStructure( NJets, Jetrate, variation.getTag() );
  return;
}

