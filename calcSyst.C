
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Analysis.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TMath.h"
#include <iostream>
#include <stdexcept>
#include <iomanip>
using std::cout;
using std::endl;
#endif

class AnalysisObject {
public:
  AnalysisObject( Int_t n ) : points(n), values(n), errors(n) {}
  virtual ~AnalysisObject() {}
  virtual TVectorD getPoints() { return points; }
  virtual TVectorD getValues() { return values; }
  virtual TVectorD getErrors() { return errors; }
  virtual TString getPointStr( Int_t ) = 0;
  virtual TVectorD getPointsCenter() = 0;
protected:
  TVectorD points;
  TVectorD values;
  TVectorD errors;
};

class TH1DAnalysisObject: public AnalysisObject {
public:
  TH1DAnalysisObject( TH1D* h, const TString & opt ) : 
    AnalysisObject( h->GetNbinsX()+1 ), hist(h) {
    if( opt.Contains( "n" ) ) normalise();
    Double_t integral= 0.0;
    Double_t integralError= 0.0;
    Int_t nbin= hist->GetNbinsX();
    for( Int_t i= 0; i < nbin; i++ ) {
      points[i]= hist->GetBinLowEdge( i+1 );
      values[i]= hist->GetBinContent( i+1 );
      errors[i]= hist->GetBinError( i+1 );
      Double_t binw= hist->GetBinWidth( i+1 );
      integral+= values[i]*binw;
      integralError+= TMath::Power( errors[i]*binw, 2 );
    }
    values[nbin]= integral;
    errors[nbin]= TMath::Sqrt( integralError );
    points[nbin]= hist->GetBinLowEdge( nbin+1 );
  }
  virtual ~TH1DAnalysisObject() {}
  virtual TString getPointStr( Int_t i ) {
    if( i < 0 || i >= points.GetNoElements() ) return "getPoint: error";
    else if( i == points.GetNoElements()-1 ) return "Integral:";
    else return Form( "%4.2f %4.2f", points[i], points[i+1] );
  }
  virtual TVectorD getPointsCenter() {
    Int_t nbin= points.GetNoElements()-1;
    TVectorD bincenters( nbin );
    for( Int_t i= 0; i < nbin; i++ ) {
      bincenters[i]= (points[i]+points[i+1])/2.0;
    }
    return bincenters;
  }
private:
  void normalise() {
    Double_t nentries= hist->GetEntries();
    for( Int_t i= 1; i < hist->GetNbinsX()+1; i++ ) {
      Double_t value= hist->GetBinContent( i );
      Double_t error= hist->GetBinError( i );
      Double_t binw= hist->GetBinWidth( i );
      Double_t scf= nentries*binw;
      hist->SetBinContent( i, value/scf );
      hist->SetBinError( i, error/scf );
    }
  }
  TH1D* hist;
};

class TGEAnalysisObject: public AnalysisObject {
public:
  TGEAnalysisObject( TGraphErrors* t, const TString & opt ) : 
    AnalysisObject( t->GetN() ), tge(t) {
    if( opt.Contains( "n" ) ) normalise();
    points.SetElements( tge->GetX() );
    values.SetElements( tge->GetY() );
    errors.SetElements( tge->GetEY() );
  }
  virtual ~TGEAnalysisObject() {}
  TString getPointStr( Int_t i ) {
    if( i < 0 || i >= points.GetNoElements() ) return "getPoint: error";
    else if( i == points.GetNoElements()-1 ) return "Norm:";
    else return Form( "%5.2f", points[i] );
  }
  virtual TVectorD getPointsCenter() {
    return points.GetSub( 0, points.GetNoElements()-2 );
  }
private:
  void normalise() {
    Int_t npoints= tge->GetN();
    Double_t* x= tge->GetX();
    Double_t* y= tge->GetY();
    Double_t* ey= tge->GetEY();
    Double_t scf= y[npoints-1];
    for( Int_t i= 0; i < npoints; i++ ) {
      tge->SetPoint( i, x[i], y[i]/scf );
      tge->SetPointError( i, 0.0, ey[i]/scf );
    }
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
					   const Analysis& anal,
					   TString opt="n" ) {
  if( opt.Contains( "v" ) ) {
    cout << "getAnalysisObjectFromFile: create " << obs << " "    
	 << anal.getTag() << endl;
  }
  TObject* tobj= getTObjFromFile( f, obs, anal );
  TString tobjclass= tobj->ClassName();
  AnalysisObject* ao= 0;
  if( tobjclass == "TH1D" ) {
    TH1D* hist= dynamic_cast<TH1D*>( tobj );
    if( hist ) ao= new TH1DAnalysisObject( hist, opt );
  }
  else if( tobjclass == "TGraphErrors" ) {
    TGraphErrors* tge= dynamic_cast<TGraphErrors*>( tobj );
    if( tge ) ao= new TGEAnalysisObject( tge, opt );
  }
  else {
    tobjclass.Prepend( "getAnalysisObjectFromFile: wrong class: " );
    throw std::logic_error( tobjclass.Data() );
  }
  return ao;
}

void insert( TVectorD& lhs, const TVectorD& rhs ) {
  lhs.ResizeTo( rhs );
  lhs= rhs;
}

struct PlotOptions {
  PlotOptions( Double_t xlo, Double_t xhi, Double_t ylo, Double_t yhi,
	       Int_t ms=20, Int_t msize=1.0 ) :
    xmin(xlo), xmax(xhi), ymin(ylo), ymax(yhi), 
    markerStyle(ms), markerSize(msize) {}
  Double_t xmin;
  Double_t xmax;
  Double_t ymin;
  Double_t ymax;
  Int_t markerStyle;
  Int_t markerSize;
};

class AnalysisObservable {
protected:
  TString obs;
  AnalysisObject* aostand;
  TVectorD points;
  TVectorD values;
  TVectorD sterrs;
  TVectorD syerrs;
public:
  AnalysisObservable( const TString& name ) : obs(name) {}
  virtual void printResults( TString opt="?",
			     Int_t width=6,
			     Int_t precision=3 ) {
    cout << "Results for " << obs << endl;
    Int_t imax= points.GetNoElements()-1;
    if( opt.Contains( "n" ) ) imax+=1;
    for( Int_t i= 0; i < imax; i++ ) {
      TString point= aostand->getPointStr(i);
      cout << point << fixed << setprecision(precision) << " ";
      cout << setw(width) << values[i] << " ";
      cout << setw(width) << sterrs[i] << " ";
      cout << setw(width) << syerrs[i] << endl;
    }
  }
  virtual void plot( const PlotOptions& po, TString opt="?" ) {
    TVectorD vx= aostand->getPointsCenter();
    TVectorD vex( vx.GetNoElements() );
    TGraphErrors* tgest= new TGraphErrors( vx, values, vex, sterrs );
    TGraphErrors* tgesy= new TGraphErrors( vx, values, vex, syerrs );
    if( opt.Contains( "s" ) ) {
      tgesy->Draw( "psame" );
    }
    else {
      tgesy->SetTitle( obs );
      tgesy->SetMinimum( po.ymin );
      tgesy->SetMaximum( po.ymax );
      tgesy->GetXaxis()->SetLimits( po.xmin, po.xmax);
      tgesy->Draw( "ap" );
    }
    tgesy->SetMarkerStyle( po.markerStyle );
    tgesy->SetMarkerSize( po.markerSize );
    tgest->Draw( "psame" );
  }
};

class LEP1AnalysisObservable : public AnalysisObservable {
public:
  LEP1AnalysisObservable( const TString& obs, TFile* f ) : AnalysisObservable(obs) {
    Analysis standardAnalysis( "data", "mt", "stand", "none", "none", "py", "bbb" );
    Analysis tcAnalysis( "data", "tc", "stand", "none", "none", "py", "bbb" );
    Analysis costt07Analysis( "data", "mt", "costt07", "none", "none", "py", "bbb" );
    Analysis nch7Analysis( "data", "mt", "nch7", "none", "none", "py", "bbb" );
    Analysis hwAnalysis( "data", "mt", "stand", "none", "none", "hw", "bbb" );
    // Analysis stand( "py", "mt", "stand", "none", "none", "py", "bbb" );
    // Analysis tc( "py", "tc", "stand", "none", "none", "py", "bbb" );
    // Analysis costt07( "py", "mt", "costt07", "none", "none", "py", "bbb" );
    // Analysis nch7( "py", "mt", "nch7", "none", "none", "py", "bbb" );
    // Analysis hw( "hw", "mt", "stand", "none", "none", "hw", "bbb" );
    try {
      aostand= getAnalysisObjectFromFile( f, obs, standardAnalysis );
      insert( points, aostand->getPoints() );
      insert( values, aostand->getValues() );
      insert( sterrs, aostand->getErrors() );
      TVectorD tctvd= getAnalysisObjectFromFile( f, obs, tcAnalysis )->getValues();
      TVectorD costt07tvd= getAnalysisObjectFromFile( f, obs, costt07Analysis )->getValues();
      TVectorD nch7tvd= getAnalysisObjectFromFile( f, obs, nch7Analysis )->getValues();
      TVectorD hwtvd= getAnalysisObjectFromFile( f, obs, hwAnalysis )->getValues();      
      TVectorD tcdelta= tctvd - values;
      TVectorD costt07delta= costt07tvd - values;
      TVectorD nch7delta= nch7tvd - values;
      TVectorD hwdelta= hwtvd - values;
      syerrs.ResizeTo( points );
      syerrs= tcdelta.Sqr() + costt07delta.Sqr() + nch7delta.Sqr() + hwdelta.Sqr();
      syerrs.Sqrt();
    }
    catch( const std::logic_error& e ) {
      cout << "Cought exception: " << e.what() << endl;
    }
  }

  // Try matrix object which throws error
  // try {
  //   Analysis matrix( "py", "mt", "stand", "nonrad", "hadron" );
  //   AnalysisObject* ao= getAnalysisObjectFromFile( f, obs, matrix );
  // }
  // catch( const std::exception& e ) {
  //   cout << "Cought exception: " << e.what() << endl;
  // }

};

void printResult( TString name="thrust", TString opt="?", Int_t width=6, Int_t precision=3 ) {
  //gROOT->LoadMacro("libNtupleReaderDict.so");
  TFile* f= new TFile( "LEP1Analysis.root" );
  LEP1AnalysisObservable obs( name, f );
  obs.printResults( opt, width, precision );
}

void plotResult( PlotOptions po, TString name="thrust", TString opt="?" ) {
  TFile* f= new TFile( "LEP1Analysis.root" );
  LEP1AnalysisObservable obs( name, f );
  obs.plot( po, opt );
}

void printResults( TString opt="?" ) {
  vector<TString> vobs;
  vobs.push_back( "thrust" );
  vobs.push_back( "antikteminR3" );
  vobs.push_back( "antiktRR3" );
  vobs.push_back( "sisconeeminR3" );
  vobs.push_back( "sisconeRR3" );
  vobs.push_back( "pxconeeminR3" );
  vobs.push_back( "pxconeRR3" );
  for( size_t iobs=0; iobs < vobs.size(); iobs++ ) {
    printResult( vobs[iobs], opt );
  }
}
