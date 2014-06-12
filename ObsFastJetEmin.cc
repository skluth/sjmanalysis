
#include "ObsFastJetEmin.hh"
#include "JetrateDataStructure.hh"
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
  vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( variation.getReco() );
  size_t n= EminFractions.size();
  vector<Double_t> NJets( n );
  TFastJet tfjakt( vtlv, "eeantikt", Rvalue );
  Double_t Evis= tfjakt.Evis();
  for( size_t i= 0; i < n; i++ ) {
    vector<TLorentzVector> incljets= tfjakt.inclusive_jets( EminFractions[i]*Evis );
    NJets[i]= incljets.size();
  }

  map<string,DataStructure*>::iterator iter= datastructures.find( variation.getTag() );
  if( iter != datastructures.end() ) {
    JetrateDataStructure* jrds= dynamic_cast<JetrateDataStructure*>( iter->second );
    if( jrds ) {
      jrds->fill( NJets, Jetrate );
    }
    else {
      std::cout << "ObsFastJetEmin::fill: dynamic_cast to JetrateDataStructure failed" << std::endl;
    }
  }
  else {
    std::cout << "ObsFastEmin::fill: " << variation.getTag() << " not found" << std::endl;
  }

}
