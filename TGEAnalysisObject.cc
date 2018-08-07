
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

TString TGEAnalysisObject::getPointStr( Int_t i, Int_t width, Int_t prec ) {
  if( i < 0 || i >= points.GetNoElements() ) return "getPoint: error";
  else {
    std::string strwidth( std::to_string( width ) );
    std::string strprec( std::to_string( prec ) );
    std::string formstr( "%"+strwidth+"."+strprec+"f" );
    return Form( formstr.c_str(), points[i] );
  }
}

TVectorD TGEAnalysisObject::getPointsCenter() {
  return points;
}

TString TGEAnalysisObject::getPointLabel( Int_t width ) {
  std::string strwidth( std::to_string( width ) );
  std::string formstr= "%-"+strwidth+"s";
  return Form( formstr.c_str(), "point" );
}
