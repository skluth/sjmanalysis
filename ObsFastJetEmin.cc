
#include "ObsFastJetEmin.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
#include <iostream>
using std::cout;
using std::endl;

ObsFastJetEmin::ObsFastJetEmin( const string& name, const string& algo,
				Double_t rval, 
				const vector<Double_t>& eminfractions, 
				const vector<Analysis>& variations ) :
  ObsJetrate( name, eminfractions ), Algorithm(algo), Rvalue(rval) {
  addAnalyses( variations );
  cout << "ObsFastJetEmin::ObsFastJetEmin: create " << getName() 
       << " with algorithm " << algo << " and R=" << rval << endl;
  printPoints();
}

void ObsFastJetEmin::fill( NtupleReader* ntr, const Analysis& variation ) {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( variation.getReco() );
  size_t n= points.size();
  vector<Double_t> NJets( n );
  TFastJet tfj( vtlv, Algorithm.c_str(), Rvalue );
  Double_t Evis= tfj.Evis();
  for( size_t i= 0; i < n; i++ ) {
    vector<TLorentzVector> incljets= tfj.inclusive_jets( points[i]*Evis );
    NJets[i]= incljets.size();
  }
  getAndFillJetrateDataStructures( NJets, variation.getTag() );
  return;
}

