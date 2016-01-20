
#include "ObsFastJetR.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
#include <iostream>
using std::cout;
using std::endl;

ObsFastJetR::ObsFastJetR( const string& name, const string& algo,
			  Double_t eminfrac,
			  const vector<Double_t>& rvalues, 
			  const vector<Analysis>& variations ) :
  ObsJetrate( name, rvalues ), Algorithm(algo), EminFraction(eminfrac) {
  addAnalyses( variations );
  cout << "ObsFastJetR::ObsFastJetR: create " << getName() 
       << " with algorithm " << algo << " and Emin/Evis=" << eminfrac << endl;
  printVectorD( "R points:", rvalues );
}

void ObsFastJetR::fill( NtupleReader* ntr, const Analysis& variation ) {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( variation.getReco() );
  size_t n= points.size();
  vector<Double_t> NJets( n );
  for( size_t i= 0; i < n; i++ ) {
    TFastJet tfj( vtlv, Algorithm.c_str(), points[i] );
    Double_t Evis= tfj.Evis();
    vector<TLorentzVector> incljets= tfj.inclusive_jets( EminFraction*Evis );
    NJets[i]= incljets.size();
  }
  getAndFillJetrateDataStructures( NJets, variation.getTag() );
  return;
}

