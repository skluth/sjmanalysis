
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "DataStructure.hh"
#include "Analysis.hh"
#include "Observable.hh"
#include "FilledObservable.hh"
#include "ObservableFactory.hh"
#include "Unfolder.hh"
#include "OutputWriter.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <string>
using std::string;
#endif

void printy23( Int_t ievnt=100, const char* filename="da91_96_200.root" ) {
  // Load libs in root before loading this macro
  // gROOT->LoadMacro("libNtupleReaderDict.so");
  NtupleReader* ntr= new NtupleReader( filename );
  if( ntr->GetEvent( ievnt ) > 0 ) {
    vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( "mt" );
    TFastJet tfjd( vtlv, "eekt" );
    TFastJet tfjj( vtlv, "jade" );
    cout << "y23 eekt values " << tfjd.ymerge( 2 ) << " " << ntr->getYmergeD( "mt", 2 ) << endl;
    cout << "y23 jade values " << tfjj.ymerge( 2 ) << " " << ntr->getYmergeE( "mt", 2 ) << endl;
  }
  else {
    cout << "event " << ievnt << " not found" << endl;
  }
  return;
}

void JetratePlot( FilledObservable* obs, const Analysis& anal ) {
  map<string,DataStructure*> data= obs->getData();
  string tag= anal.getTag();
  DataStructure* jrds= data[tag];
  string txt= obs->getName() + " " + tag;
  vector<Double_t> points= jrds->getPoints();
  vector<Double_t> values= jrds->getValues();
  vector<Double_t> errors= jrds->getErrors();
  Int_t n= points.size();
  Double_t* xerrors= new Double_t[n];
  for( Int_t i= 0; i < n; i++ ) {
    xerrors[i]= 0.0;
  }
  TGraphErrors* tg= new TGraphErrors( n, &(points[0]), &(values[0]), xerrors, &(errors[0]) );
  tg->SetMarkerStyle( 20 );
  tg->SetMarkerSize( 1.0 );
  tg->SetTitle( txt.c_str() );
  tg->Draw( "ap" );
  return;
}

void processDurhamJade( const vector<TLorentzVector>& vtlv, Double_t& y23d, Double_t& y23j ) {
  TFastJet tfjDurham( vtlv, "eekt" );
  y23d= tfjDurham.ymerge(2);
  TFastJet tfjJade( vtlv, "jade" );
  y23j= tfjJade.ymerge(2);
  return;
}

vector<FilledObservable*> getFilled( const vector<Observable*>& vobs ) {
  vector<FilledObservable*> vfobs;
  for( size_t iobs= 0; iobs < vobs.size(); iobs++ ) {
    vector<FilledObservable*> vfobspart= vobs[iobs]->getFilledObservables();
    vfobs.insert( vfobs.end(), vfobspart.begin(), vfobspart.end() );
  }
  return vfobs;
}

FilledObservable* findFilledObservable( TString name, 
					const vector<FilledObservable*>& vfobs ) {
  FilledObservable* fobs= 0;
  for( size_t ifobs= 0; ifobs < vfobs.size(); ifobs++ ) {
    TString obsname= vfobs[ifobs]->getName();
    if( obsname == name ) {
      fobs= vfobs[ifobs];
      break;
    }
  }
  if( fobs == 0 ) cout << "FilledObservable " << name << " not found" << endl;
  return fobs;
}

void makeplots( Int_t maxevt=1000, 
		const char* pyfilename="mc5025_1_200.root" ) {

  // Load libs in root before loading this macro
  // gROOT->LoadMacro("libNtupleReaderDict.so");
  
  TH1F* y23dmt= new TH1F( "y23dmt", "y23 D det mt", 10, 0.0, 5.0 );
  TH1F* y23dmtnr= new TH1F( "y23dmtnr", "y23 D det mt non-rad", 10, 0.0, 5.0 );
  TH1F* y23dhad= new TH1F( "y23dhad", "y23 D had", 10, 0.0, 5.0 );
  TH1F* y23dhadsel= new TH1F( "y23dhadsel", "y23 D had sel", 10, 0.0, 5.0 );
  TH2F* y23dmthad= new TH2F( "y23dmthad", "y23 D mt vs had non-rad", 10, 0.0, 5.0, 10, 0.0, 5.0 );

  TH1F* y23jmt= new TH1F( "y23jmt", "y23 J det mt", 10, 0.0, 5.0 );
  TH1F* y23jmtnr= new TH1F( "y23jmtnr", "y23 J det mt non-rad", 10, 0.0, 5.0 );
  TH1F* y23jhad= new TH1F( "y23jhad", "y23 J had", 10, 0.0, 5.0 );
  TH1F* y23jhadsel= new TH1F( "y23jhadsel", "y23 J had sel", 10, 0.0, 5.0 );
  TH2F* y23jmthad= new TH2F( "y23jmthad", "y23 J mt vs had non-rad", 10, 0.0, 5.0, 10, 0.0, 5.0 );

  vector<Analysis> allAnalyses;
  allAnalyses.push_back( Analysis( "py", "mt", "stand" ) );
  allAnalyses.push_back( Analysis( "py", "hadron", "none", "nonrad" ) );

  vector<string> obsnames;
  obsnames.push_back( "antiktemin" );
  obsnames.push_back( "antiktR" );
  ObservableFactory obsfac;
  vector<Observable*> vobs= obsfac.createObservables( obsnames, allAnalyses );

  vector<string> selections;
  selections.push_back( "stand" );
  selections.push_back( "costt07" );
  selections.push_back( "nch7" );
  selections.push_back( "nonrad" );
  selections.push_back( "both" );
  map<TString,int> Nselected;
  for( size_t isel= 0; isel < selections.size(); isel++ ) {
    string sel= selections[isel];
    Nselected[sel]= 0;
  }

  NtupleReader* ntr= new NtupleReader( pyfilename );
  Int_t nevnt= ntr->GetNumberEntries();
  for( Int_t ievnt= 0; ievnt < TMath::Min( nevnt, maxevt ); ievnt++ ) {

    if( ntr->GetEvent( ievnt ) == 0 ) {
      cout << "event " << ievnt << " not found, skip" << endl;
      continue;
    }

    // Fill Observables:
    for( size_t i= 0; i < allAnalyses.size(); i++ ) {
      string cuts= allAnalyses[i].getCuts();
      string mccuts= allAnalyses[i].getMccuts();
      map<string,Bool_t> LEP1selections= ntr->LEP1Selections();
      bool MCnonrad= ntr->MCNonRad();
      if( ( cuts == "none" || LEP1selections[cuts] ) &&
	  ( mccuts == "none" || MCnonrad ) ) {
	for( size_t j= 0; j < vobs.size(); j++ ) {
	  vobs[j]->fill( ntr, allAnalyses[i] );
	}
      }
    }

    Double_t y23dmtfj= -1.0, y23jmtfj= -1.0;
    if( ntr->LEP1Preselection() ) {

      map<string,bool> isSelected= ntr->LEP1Selections();
      if( isSelected["stand"] ) {
	vector<TLorentzVector> vtlvmt= ntr->GetLorentzVectors( "mt" );
	processDurhamJade( vtlvmt, y23dmtfj, y23jmtfj );
	y23dmt->Fill( -TMath::Log10( y23dmtfj ) );
	y23jmt->Fill( -TMath::Log10( y23jmtfj ) );

	Nselected["stand"]++;
	if( isSelected["costt07"] ) Nselected["costt07"]++;
	if( isSelected["nch7"] ) Nselected["nch7"]++;

      // Float_t y23dold= ntr->dmt_ymerge(2);
      // Float_t y23jold= ntr->emt_ymerge(2);
      // cout << y23dold << " " << y23dmtfj << " "
      // 	   << y23jold << " " << y23jmtfj << endl;
      // Get jet 4-vectors:
      // vector<TLorentzVector> jets= tfjmt.exclusive_jets( 4 );
      // cout << "excl. jets " << ievnt<< " " << jets.size() << endl;
      // for( Int_t i= 0; i < 4; i++ ) {
      // 	jets[i].Print();
      // }
      }

    }
    Double_t y23dhadfj= -1.0, y23jhadfj= -1.0;
    if( ntr->MCNonRad() ) {
      Nselected["nonrad"]++;
      vector<TLorentzVector> vtlvhad= ntr->GetLorentzVectors( "hadron" );
      processDurhamJade( vtlvhad, y23dhadfj, y23jhadfj );
      y23dhad->Fill( -TMath::Log10( y23dhadfj ) );
      y23jhad->Fill( -TMath::Log10( y23jhadfj ) );
    }
    if( ntr->LEP1Selection() && ntr->MCNonRad() ) {
      Nselected["both"]++;
      y23dmthad->Fill( -TMath::Log10( y23dhadfj ), -TMath::Log10( y23dmtfj ) );
      y23dhadsel->Fill( -TMath::Log10( y23dhadfj ) );
      y23dmtnr->Fill( -TMath::Log10( y23dmtfj ) );
      y23jmthad->Fill( -TMath::Log10( y23jhadfj ), -TMath::Log10( y23jmtfj ) );
      y23jhadsel->Fill( -TMath::Log10( y23jhadfj ) );
      y23jmtnr->Fill( -TMath::Log10( y23jmtfj ) );
    }
  }
  cout << "Total number of events " << nevnt << endl;
  cout << "LEP1 selected events " << Nselected["stand"] << endl;
  cout << "LEP1 selected events costt<0.7 " << Nselected["costt07"] << endl;
  cout << "LEP1 selected events nch>=7 " << Nselected["nch7"] << endl;
  cout << "MC non-rad selected events " <<  Nselected["nonrad"] << endl;
  cout << "LEP1 and MC non-rad selected events " << Nselected["both"] << endl;

  // Get filled observables for plots:
  vector<FilledObservable*> vfobs= getFilled( vobs );
  FilledObservable* antikteminR3= findFilledObservable( "antikteminR3", vfobs );
  FilledObservable* antiktRR3= findFilledObservable( "antiktRR3", vfobs );
  if( antikteminR3 == 0 || antiktRR3 == 0 ) return;

  // Plots:
  TCanvas* canvd= new TCanvas( "canvd", "y23 Durham plots", 900, 900 );
  canvd->Divide( 2, 2 );
  canvd->cd( 1 );
  y23dmt->Draw();
  canvd->cd( 2 );
  y23dhadsel->Draw();
  canvd->cd( 3 );
  y23dmtnr->Draw();
  canvd->cd( 4 );
  y23dmthad->Draw( "text" );

  TCanvas* canvj= new TCanvas( "canvj", "y23 Jade plots", 900, 900 );
  canvj->Divide( 2, 2 );
  canvj->cd( 1 );
  y23jmt->Draw();
  canvj->cd( 2 );
  y23jhadsel->Draw();
  canvj->cd( 3 );
  y23jmtnr->Draw();
  canvj->cd( 4 );
  y23jmthad->Draw( "box" );

  Analysis detectorLevel( "py", "mt", "stand" );
  Analysis hadronLevel( "py", "hadron", "none", "nonrad" );
  TCanvas* canva= new TCanvas( "canva", "ee anti-kt plots", 900, 900 );
  canva->Divide( 2, 2 );
  canva->cd( 1 );
  JetratePlot( antikteminR3, hadronLevel );
  canva->cd( 2 );
  JetratePlot( antiktRR3, hadronLevel );
  canva->cd( 3 );
  JetratePlot( antikteminR3, detectorLevel );
  canva->cd( 4 );
  JetratePlot( antiktRR3, detectorLevel );

  // canv->Print( "plot.pdf" );
    
  return;

}

