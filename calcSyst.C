
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Analysis.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TMath.h"
#include <iostream>
#include <stdexcept>
using std::cout;
using std::endl;
#endif

class AnalysisObject {
public:
  virtual TVectorD getTVectorD( TString opt="v" ) = 0;
};

class TH1DAnalysisObject: public AnalysisObject {
public:
  TH1DAnalysisObject( TH1D* h ) : hist(h) {}
  TVectorD getTVectorD( TString opt="v" ) {
    TH1D* tmphist= hist;
    if( opt.Contains( "n" ) ) tmphist= normalise();
    Int_t nbin= tmphist->GetNbinsX();
    TVectorD tvd( nbin+1 );
    Double_t integral= 0.0;
    for( Int_t i= 0; i < nbin; i++ ) {
      Double_t binw= tmphist->GetBinWidth( i+1 );
      if( opt.Contains( "v" ) ) {
	tvd[i]= tmphist->GetBinContent( i+1 );
	integral+= tvd[i]*binw;
      }
      else if( opt.Contains( "e" ) ) {
	tvd[i]= tmphist->GetBinError( i+1 );
	integral+= TMath::Power( tvd[i]*binw, 2 );
      }
    }
    if( opt.Contains( "v" ) ) tvd[nbin]= integral;
    else if( opt.Contains( "e" ) ) tvd[nbin]= TMath::Sqrt( integral );
    if( opt.Contains( "n" ) ) delete tmphist;
    return tvd; 
  }
private:
  TH1D* normalise() {
    TH1D* nhist= new TH1D( *hist );
    Double_t nentries= hist->GetEntries();
    for( Int_t i= 1; i < hist->GetNbinsX()+1; i++ ) {
      Double_t binw= hist->GetBinWidth( i );
      Double_t scf= nentries*binw;
      Double_t value= hist->GetBinContent( i );
      Double_t error= hist->GetBinError( i );
      nhist->SetBinContent( i, value/scf );
      nhist->SetBinError( i, error/scf );
    }
    return nhist;
  }
  TH1D* hist;
};

class TGEAnalysisObject: public AnalysisObject {
public:
  TGEAnalysisObject( TGraphErrors* t ) : tge(t) {}
  TVectorD getTVectorD( TString opt="v" ) {
    TGraphErrors* tmptge= tge;
    if( opt.Contains( "n" ) ) tmptge= normalise();
    TVectorD tvd( tmptge->GetN() );
    if( opt.Contains( "v" ) ) tvd.SetElements( tmptge->GetY() );
    else if( opt.Contains( "e" ) ) tvd.SetElements( tmptge->GetEY() );
    if( opt.Contains( "n" ) ) delete tmptge;
    return tvd;
  }
private:
  TGraphErrors* normalise() {
    Int_t npoint= tge->GetN();
    TGraphErrors* ntge= new TGraphErrors( *tge );
    Double_t* abs= tge->GetX();
    Double_t* values= tge->GetY();
    Double_t* errors= tge->GetEY();
    Double_t scf= values[npoint-1];
    for( Int_t i= 0; i < npoint; i++ ) {
      ntge->SetPoint( i, abs[i], values[i]/scf );
      ntge->SetPointError( i, 0.0, errors[i]/scf );
    }
    return ntge;
  }
  TGraphErrors* tge;
};

TObject* getTObjFromFile( TFile* f, const TString& obs, const Analysis& anal ) {
  TString tag= anal.getTag();
  TString key= obs + " " + tag;
  TObject* tobj= (TObject*) f->Get( key );
  if( tobj == 0 ) {
    key.Prepend( "getTObjFromFile: no object with key: " );
    throw std::logic_error( key.Data() );
  }
  return tobj;
}

AnalysisObject* getAnalysisObjectFromFile( TFile* f, const TString& obs, 
					   const Analysis& anal ) {
  cout << "getAnalysisObjectFromFile: create " << obs << " "    
       << anal.getTag() << endl;
  TObject* tobj= getTObjFromFile( f, obs, anal );
  TString tobjclass= tobj->ClassName();
  AnalysisObject* ao= 0;
  if( tobjclass == "TH1D" ) {
    TH1D* hist= dynamic_cast<TH1D*>( tobj );
    if( hist ) ao= new TH1DAnalysisObject( hist );
  }
  else if( tobjclass == "TGraphErrors" ) {
    TGraphErrors* tge= dynamic_cast<TGraphErrors*>( tobj );
    if( tge ) ao= new TGEAnalysisObject( tge );
  }
  else {
    tobjclass.Prepend( "getAnalysisObjectFromFile: wrong class: " );
    throw std::logic_error( tobjclass.Data() );
  }
  return ao;
}


void printResults( TString obs="thrust" ) {

  //gROOT->LoadMacro("libNtupleReaderDict.so");

  Analysis stand( "data", "mt", "stand", "none", "none", "py", "bbb" );
  Analysis tc( "data", "tc", "stand", "none", "none", "py", "bbb" );
  Analysis costt07( "data", "mt", "costt07", "none", "none", "py", "bbb" );
  Analysis nch7( "data", "mt", "nch7", "none", "none", "py", "bbb" );
  Analysis hw( "data", "mt", "stand", "none", "none", "hw", "bbb" );

  // Analysis stand( "py", "mt", "stand", "none", "none", "py", "bbb" );
  // Analysis tc( "py", "tc", "stand", "none", "none", "py", "bbb" );
  // Analysis costt07( "py", "mt", "costt07", "none", "none", "py", "bbb" );
  // Analysis nch7( "py", "mt", "nch7", "none", "none", "py", "bbb" );
  // Analysis hw( "hw", "mt", "stand", "none", "none", "hw", "bbb" );

  TFile* f= new TFile( "LEP1Analysis.root" );

  try {

    TVectorD standtvd= getAnalysisObjectFromFile( f, obs, stand )->getTVectorD( "vn" );
    TVectorD standerrtvd= getAnalysisObjectFromFile( f, obs, stand )->getTVectorD( "en" );
    TVectorD tctvd= getAnalysisObjectFromFile( f, obs, tc )->getTVectorD( "vn" );
    TVectorD costt07tvd= getAnalysisObjectFromFile( f, obs, costt07 )->getTVectorD( "vn" );
    TVectorD nch7tvd= getAnalysisObjectFromFile( f, obs, nch7 )->getTVectorD( "vn" );
    TVectorD hwtvd= getAnalysisObjectFromFile( f, obs, hw )->getTVectorD( "vn" );

    TVectorD tcdelta= tctvd - standtvd;
    TVectorD costt07delta= costt07tvd - standtvd;
    TVectorD nch7delta= nch7tvd - standtvd;
    TVectorD hwdelta= hwtvd - standtvd;

    TVectorD systerr= tcdelta.Sqr() + costt07delta.Sqr() + nch7delta.Sqr() + hwdelta.Sqr();
    systerr.Sqrt();

    for( Int_t i= 0; i < systerr.GetNoElements(); i++ ) {
      cout << i << " " << standtvd[i] << " +/- " << standerrtvd[i]
	   << " +/- " << systerr[i] << endl;
    }

  }

  catch( const std::exception& e ) {
    cout << "Cought exception: " << e.what() << endl;
  }

  // Try matrix object which throuws error
  try {
    Analysis matrix( "py", "mt", "stand", "nonrad", "hadron" );
    AnalysisObject* ao= getAnalysisObjectFromFile( f, obs, matrix );
  }
  catch( const std::exception& e ) {
    cout << "Cought exception: " << e.what() << endl;
  }

  return;
  
}


