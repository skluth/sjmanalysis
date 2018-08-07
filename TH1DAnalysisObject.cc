
#include "TH1DAnalysisObject.hh"

#include "Normalisation.hh"

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

#include <iostream>

TH1DAnalysisObject::TH1DAnalysisObject( TH1D* h, TH2D* h2d ) :
  AnalysisObject( h->GetNbinsX()+1 ), hist(h), hist2d(h2d) {
  Double_t integral= 0.0;
  Double_t integralError= 0.0;
  UInt_t nbin= hist->GetNbinsX();
  for( UInt_t i= 0; i < nbin; i++ ) {
    points[i]= hist->GetBinLowEdge( i+1 );
    Double_t binw= hist->GetBinWidth( i+1 );
    if( IsNormalised( *hist ) ) {
      values[i]= hist->GetBinContent( i+1 );
      errors[i]= hist->GetBinError( i+1 );
    }
    else {
      Double_t norm= hist->GetEntries()*binw;
      values[i]= hist->GetBinContent( i+1 )/norm;
      errors[i]= hist->GetBinError( i+1 )/norm;
    }
    integral+= values[i]*binw;
    integralError+= TMath::Power( errors[i]*binw, 2 );
  }
  values[nbin]= integral;
  errors[nbin]= TMath::Sqrt( integralError );
  points[nbin]= hist->GetBinLowEdge( nbin+1 );
  if( hist2d != 0 ) {
    UInt_t ndim= hist2d->GetNbinsX();
    errorMatrix.ResizeTo( ndim, ndim );
    for( UInt_t i= 0; i < ndim; i++ ) {
      for( UInt_t j= 0; j < ndim; j++ ) {
	errorMatrix[i][j]= hist2d->GetBinContent( i+1, j+1 );
      }
    }
  }
}

TString TH1DAnalysisObject::getPointStr( Int_t i, Int_t width, Int_t prec ) {
  if( i < 0 || i >= points.GetNoElements() ) return "getPoint: error";
  else if( i == points.GetNoElements()-1 ) {
    std::string strwidth( std::to_string( 2*width+1 ) );
    std::string formstr( "%-"+strwidth+"s" );
    return Form( formstr.c_str(), "Integral:" );
  }
  else {
    std::string strwidth( std::to_string( width ) );
    std::string strprec( std::to_string( prec ) );
    std::string formstr( "%"+strwidth+"."+strprec+"f %"+strwidth+"."+strprec+"f" );
    return Form( formstr.c_str(), points[i], points[i+1] );
  }
}

TVectorD TH1DAnalysisObject::getPointsCenter() {
  Int_t nbin= points.GetNoElements()-1;
  TVectorD bincenters( nbin );
  for( Int_t i= 0; i < nbin; i++ ) {
    bincenters[i]= (points[i]+points[i+1])/2.0;
  }
  return bincenters;
}

TString TH1DAnalysisObject::getPointLabel( Int_t width ) {
  std::string strwidth( std::to_string( width ) );
  std::string formstr( "%-"+strwidth+"s %-"+strwidth+"s" );
  return Form( formstr.c_str(), "lo", "hi" );
}

TVectorD TH1DAnalysisObject::getErrors( const TString & opt ) {
  TVectorD result= AnalysisObject::getErrors( opt );
  // With full stat. error matrix error on integral is 0, since the
  // integral is constrained by the normalisation
  if( opt.Index( "m" ) >= 0 ) {
    if( errorMatrix.GetNoElements() > 0 ) {
      result[result.GetNrows()-1]= 0.0;
    }
    else {
      std::cout << "TH1DAnalysisObject::getErrors: error matrix empty" << std::endl;
    }
  }
  return result;
}
