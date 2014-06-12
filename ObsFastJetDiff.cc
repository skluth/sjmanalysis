
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
  map<string,DataStructure*>::iterator iter= datastructures.find( variation.getTag() );
  if( iter != datastructures.end() ) {
    DifferentialDataStructure* dds= 
      dynamic_cast<DifferentialDataStructure*>( iter->second );
    if( dds ) {
      dds->fill( -TMath::Log10( yflip ) );
    }
    else {
      std::cout << "ObsFastJetDiff::fill: dynamic_cast to DifferentialDataStructure failed" 
		<< std::endl;
    }
  }
  else {
    std::cout << "ObsFastJetDiff::fill: analysis " << variation.getTag() << " not found" 
	      << std::endl;
  }
  return;
}

