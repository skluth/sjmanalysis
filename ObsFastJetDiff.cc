
#include "ObsFastJetDiff.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
#include "DataStructure.hh"
#include "DifferentialDataStructure.hh"

#include <iostream>
using std::cout;
using std::endl;
#include "TMath.h"

ObsFastJetDiff::ObsFastJetDiff( const string& name, const string& algo, Int_t njet, 
				const vector<Double_t>& bins, 
				const vector<Analysis>& variations ) :
  ObsDifferential(name,bins,variations), Algorithm(algo), Jetlower(njet) {}

void ObsFastJetDiff::fill( NtupleReader* ntr, const Analysis& variation ) {
  vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( variation.getReco() );
  TFastJet tfj( vtlv, Algorithm.c_str() );
  Double_t yflip= tfj.ymerge( Jetlower );
  getAndFillDifferentialDataStructure( -TMath::Log10( yflip ), variation.getTag() );
  return;
}

