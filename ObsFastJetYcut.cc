
#include "ObsFastJetYcut.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
#include "TMath.h"
#include <iostream>
using std::cout;
using std::endl;

ObsFastJetYcut::ObsFastJetYcut( const string& name, const string& algo,
				const vector<Double_t>& ycutpoints, 
				const vector<Analysis>& variations ) :
  ObsJetrate( name, ycutpoints ), Algorithm( algo ) {
  addAnalyses( variations );
  cout << "ObsFastJetYcut::ObsFastJetYcut: create " << getName() 
       << " with algorithm " << algo << endl;
  printPoints();
}

void ObsFastJetYcut::fill( NtupleReader* ntr, const Analysis& variation ) {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( variation.getReco() );
  size_t n= points.size();
  vector<Double_t> NJets( n );
  TFastJet tfj( vtlv, Algorithm.c_str() );
  for( size_t i= 0; i < n; i++ ) {
    NJets[i]= tfj.njets( TMath::Power( 10.0, -points[i] ) );
  }
  getAndFillJetrateDataStructures( NJets, variation.getTag() );
  return;
}

