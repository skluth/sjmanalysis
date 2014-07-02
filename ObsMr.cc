
#include "ObsMr.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
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
  printBinedges();
}

void ObsMr::fill( NtupleReader* ntr, const Analysis& variation ) {
  string reco= variation.getReco();
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( reco );
  TFastJet tfj( vtlv, "eekt" );
  vector<TLorentzVector> jets= tfj.exclusive_jets( 2 );
  if( tfj.ymerge( 3 ) > y34cut and 
      tfj.ymerge( 3 )/tfj.ymerge( 2 ) > y34y23cut ) {
    Double_t mrvalue= TMath::Min( jets[1].M2(), jets[0].M2() ) /
      TMath::Max( jets[1].M2(), jets[0].M2() );
    string tag= variation.getTag();
    getAndFillDifferentialDataStructure( mrvalue, tag, datastructures );
  }
}

