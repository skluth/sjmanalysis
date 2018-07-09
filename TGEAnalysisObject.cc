
#include "TGEAnalysisObject.hh"

#include "Normalisation.hh"

#include "TGraphErrors.h"

#include <iostream>

TGEAnalysisObject::TGEAnalysisObject( TGraphErrors* t ) : 
  AnalysisObject( IsNormalised( *t ) ? t->GetN() : t->GetN()-1 ), tge(t) {
  points.SetElements( tge->GetX() );
  if( IsNormalised( *tge ) ) {
    values.SetElements( tge->GetY() );
    errors.SetElements( tge->GetEY() );
  }
  else {
    Int_t npoints= tge->GetN()-1;
    Double_t* tgevalues= tge->GetY();
    Double_t* tgeerrors= tge->GetEY();
    Double_t nentries= tgevalues[npoints];
    for( Int_t i= 0; i < npoints; i++ ) {
      values[i]= tgevalues[i]/nentries;
      errors[i]= tgeerrors[i]/nentries;
    }
  }
}

TString TGEAnalysisObject::getPointStr( Int_t i ) {
  if( i < 0 || i >= points.GetNoElements() ) return "getPoint: error";
  else return Form( "%5.2f", points[i] );
}

TVectorD TGEAnalysisObject::getPointsCenter() {
  return points;
}

TString TGEAnalysisObject::getPointLabel() {
  return Form( "%-5s", "point" );
}
