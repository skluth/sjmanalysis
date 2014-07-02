
#include "ObsMr.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
#include "MatrixDataStructure.hh"
#include "TMath.h"
#include <iostream>
using std::cout;
using std::endl;

ObsMr::ObsMr( const vector<Double_t>& bins, 
	      const vector<Analysis>& variations,
	      Double_t y34c, Double_t y34y23c ) :
  ObsDifferential( "mr", bins ), y34cut(y34c), y34y23cut(y34y23c) {
  addAnalyses( variations );
  cout << "ObsMr::ObsMr: create " << getName() << " with cuts y34cut=" 
       << y34cut << " and y34y23cut=" << y34y23cut << endl;
  printVectorD( "Binedges:", bins );
}

void ObsMr::fill( NtupleReader* ntr, const Analysis& variation ) {
  string reco= variation.getReco();
  Double_t mrvalue= getMrvalue( ntr, reco );
  string tag= variation.getTag();
  getAndFillDifferentialDataStructure( mrvalue, tag, datastructures );
  string reco2= variation.getReco2();
  if( reco2 != "none" and ntr->isMC() ) {
    Double_t MCmrvalue= getMrvalue( ntr, reco2 );
    matrices[tag]->fill( MCmrvalue, mrvalue );
  }
}

Double_t ObsMr::getMrvalue( NtupleReader* ntr, const string& reco ) {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( reco );
  TFastJet tfj( vtlv, "eekt" );
  Double_t mrvalue= -1.0;
  if( tfj.ymerge( 3 ) > y34cut and 
      tfj.ymerge( 3 )/tfj.ymerge( 2 ) > y34y23cut ) {
    vector<TLorentzVector> jets= tfj.exclusive_jets( 2 );
    mrvalue= TMath::Min( jets[1].M2(), jets[0].M2() ) /
      TMath::Max( jets[1].M2(), jets[0].M2() );
  }
  return mrvalue;
}
