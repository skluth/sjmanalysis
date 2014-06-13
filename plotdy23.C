
#if defined(__CINT__)
//#include "Observable.hh"
#include "ObsThrust.hh"
//#include "ObsFastJetR.hh"
//#include "ObsFastJetEmin.hh"
//#include "ObsFastJetDiff.hh"
//#include <vector>
//using std::vector;
#endif

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TLorentzVector.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

#include "JetrateDataStructure.hh"
#include "Analysis.hh"
#include "Observable.hh"
#include "ObsThrust.hh"
#include "ObsFastJetR.hh"
#include "ObsFastJetEmin.hh"
#include "ObsFastJetDiff.hh"
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


void printy23( Int_t ievnt=100, const char* filename="da130_95_200.root" ) {

  // Load libs in root before loading this macro
  // gROOT->LoadMacro("libNtupleReaderDict.so");

  NtupleReader* ntr= new NtupleReader( filename );
  if( ntr->GetEvent( ievnt ) > 0 ) {
    vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( "mt" );
    TFastJet tfj( vtlv, "eekt" );
    cout << "y23 values " << tfj.ymerge(2) << " " << ntr->dmt_ymerge(2) << endl;
  }
  else {
    cout << "event " << ievnt << " not found" << endl;
  }
  return;

}


void JetratePlot( Observable* obs, const Analysis& anal ) {
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

void processAnalyses( const vector<Analysis>& analyses,
		      const vector<Observable*>& vobs,
		      const string& filename,
		      Int_t maxevt ) {
  cout << "processAnalyses: file " << filename << ", analyses:" << endl;
  for( size_t i= 0; i < analyses.size(); i++ ) {
    cout << analyses[i].getTag() << endl;
  }
  NtupleReader* ntr= new NtupleReader( filename.c_str() );
  Int_t nevnt= ntr->GetNumberEntries();
  for( Int_t ievnt= 0; ievnt < TMath::Min( nevnt, maxevt ); ievnt++ ) {
    if( ntr->GetEvent( ievnt ) == 0 ) {
      cout << "event " << ievnt << " not found, skip" << endl;
      continue;
    }
    map<string,Bool_t> selections= ntr->LEP1Selections();
    bool MCnonrad= ntr->MCNonRad();
    for( size_t i= 0; i < analyses.size(); i++ ) {
      string cuts= analyses[i].getCuts();
      string mccuts= analyses[i].getMccuts();
      if( ( cuts == "none" || selections[cuts] ) &&
	  ( mccuts == "none" || MCnonrad ) ) {
	for( size_t j= 0; j < vobs.size(); j++ ) {
	  vobs[j]->fill( ntr, analyses[i] );
	}
      }
    }
  }
  delete ntr;
  return;
}

void processUnfolding( const vector<Analysis>& measuredAnalyses, string unfoldsource,
		       vector<Observable*>& vobs ) {
  cout << "processUnfolding: bin-by-bin unfolding for analyses:" << endl;
  Analysis hadronlevel( unfoldsource, "hadron", "none", "nonrad" );
  for( size_t ianal= 0; ianal < measuredAnalyses.size(); ianal++ ) {
    Analysis measured= measuredAnalyses[ianal];
    cout << measured.getTag() << endl;
    Analysis measuredMC( unfoldsource, measured.getReco(), measured.getCuts() );
    Unfolder unfolder( measured, measuredMC, hadronlevel );
    for( size_t iobs= 0; iobs < vobs.size(); iobs++ ) {
      unfolder.unfold( vobs[iobs] );
    }
  }
  cout << "Hadron level: " << hadronlevel.getTag() << endl;
  return;
}

void LEP1Analysis( Int_t maxevt=1000, 
		   const char* datafilename="da91_96_200.root",
		   const char* pyfilename="mc5025_1_200.root", 
		   const char* hwfilename="mc12406_1_200.root" ) {

  // Load libs in root before loading this macro
  // gROOT->LoadMacro("libNtupleReaderDict.so");

  // Define analysis variations:
  vector<Analysis> measuredAnalyses;
  measuredAnalyses.push_back( Analysis( "data", "mt", "stand" ) );
  measuredAnalyses.push_back( Analysis( "data", "tc", "stand" ) );
  measuredAnalyses.push_back( Analysis( "data", "mt", "costt07" ) );
  measuredAnalyses.push_back( Analysis( "data", "mt", "nch7" ) );
  vector<Analysis> pyAnalyses;
  pyAnalyses.push_back( Analysis( "py", "mt", "stand" ) );
  pyAnalyses.push_back( Analysis( "py", "tc", "stand" ) );
  pyAnalyses.push_back( Analysis( "py", "mt", "costt07" ) );
  pyAnalyses.push_back( Analysis( "py", "mt", "nch7" ) );
  pyAnalyses.push_back( Analysis( "py", "hadron", "none", "nonrad" ) );
  pyAnalyses.push_back( Analysis( "py", "mt", "stand", "nonrad" ) );
  vector<Analysis> hwAnalyses;
  hwAnalyses.push_back( Analysis( "hw", "mt", "stand" ) );
  hwAnalyses.push_back( Analysis( "hw", "hadron", "none", "nonrad" ) );
  vector<Analysis> allAnalyses( measuredAnalyses );
  allAnalyses.insert( allAnalyses.end(), pyAnalyses.begin(), pyAnalyses.end() );
  allAnalyses.insert( allAnalyses.end(), hwAnalyses.begin(), hwAnalyses.end() );

  // Define observables:
  // Thrust
  vector<Double_t> thrustbins( 13 );
  thrustbins[0]= 0.00;
  thrustbins[1]= 0.01;
  thrustbins[2]= 0.02;
  thrustbins[3]= 0.03;
  thrustbins[4]= 0.04;
  thrustbins[5]= 0.05;
  thrustbins[6]= 0.07;
  thrustbins[7]= 0.09;
  thrustbins[8]= 0.12;
  thrustbins[9]= 0.15;
  thrustbins[10]= 0.22;
  thrustbins[11]= 0.30;
  thrustbins[12]= 0.50;
  ObsThrust* obsthrust= new ObsThrust( thrustbins, allAnalyses );
  // y23 Durham and JadeE0
  vector<Double_t> y23bins( 11 );
  for( size_t i= 0; i < y23bins.size(); i++ ) {
    y23bins[i]= 0.5*i;
  }
  ObsFastJetDiff* obsy23d= new ObsFastJetDiff( "y23d", "eekt", 2, y23bins, allAnalyses );
  ObsFastJetDiff* obsy23j= new ObsFastJetDiff( "y23j", "jade", 2, y23bins, allAnalyses );
  // anti-kt 3-jet rate vs Emin/Evis for R=0.7
  vector<Double_t> eminFraction( 9 );
  eminFraction[0]= 0.02;
  eminFraction[1]= 0.04;
  eminFraction[2]= 0.06;
  eminFraction[3]= 0.08;
  eminFraction[4]= 0.10;
  eminFraction[5]= 0.12;
  eminFraction[6]= 0.14;
  eminFraction[7]= 0.16;
  eminFraction[8]= 0.18;
  ObsFastJetEmin* obsakteminR3= new ObsFastJetEmin( "antikteminR3", "eeantikt", 3, 0.7, 
						    eminFraction, allAnalyses );
  // anti-kt 3-jet rate vs R for Emin/Evis= 0.06
  vector<Double_t> Rvalues( 8 );
  Rvalues[0]= 0.2;
  Rvalues[1]= 0.4;
  Rvalues[2]= 0.6;
  Rvalues[3]= 0.7;
  Rvalues[4]= 0.8;
  Rvalues[5]= 1.0;
  Rvalues[6]= 1.2;
  Rvalues[7]= 1.4;
  ObsFastJetR* obsaktRR3= new ObsFastJetR( "antiktRR3", "eeantikt", 3, 0.06, 
					   Rvalues, allAnalyses );
  // All observables in vector
  vector<Observable*> vobs;
  vobs.push_back( obsthrust );
  vobs.push_back( obsy23d );
  vobs.push_back( obsy23j );
  vobs.push_back( obsakteminR3 );
  vobs.push_back( obsaktRR3 );
  cout << "LEP1Analysis: observables:" << endl;
  for( size_t iobs= 0; iobs < vobs.size(); iobs++ ) {
    cout << vobs[iobs]->getName() << " ";
  }
  cout << endl;

  // Fill from data and mc (PYTHIA and HERWIG) ntuples:
  processAnalyses( measuredAnalyses, vobs, datafilename, maxevt );
  processAnalyses( pyAnalyses, vobs, pyfilename, maxevt );
  processAnalyses( hwAnalyses, vobs, hwfilename, maxevt );

  // Unfolding bin-by-bin
  // PYHTHIA based
  processUnfolding( measuredAnalyses, "py", vobs );
  // HERWIG based for systematic
  vector<Analysis> measuredAnalysesHw;
  measuredAnalysesHw.push_back( measuredAnalyses[0] );
  processUnfolding( measuredAnalysesHw, "hw", vobs );
  // MC detector level with MC as cross check for PYTHIA and HERWIG
  vector<Analysis> measuredPyAnalyses;
  measuredPyAnalyses.push_back( Analysis( "py", "mt", "stand" ) );
  measuredPyAnalyses.push_back( Analysis( "py", "tc", "stand" ) );
  measuredPyAnalyses.push_back( Analysis( "py", "mt", "costt07" ) );
  measuredPyAnalyses.push_back( Analysis( "py", "mt", "nch7" ) );
  processUnfolding( measuredPyAnalyses, "py", vobs );
  vector<Analysis> measuredHwAnalyses;
  measuredHwAnalyses.push_back( Analysis( "hw", "mt", "stand" ) );
  processUnfolding( measuredHwAnalyses, "hw", vobs );

  // Normalise and calculate stat errors, print
  for( size_t i= 0; i < vobs.size(); i++ ) {
    vobs[i]->finalise();
    vobs[i]->print();
  }

  // no systematics yet

  OutputWriter writer( "LEP1Analysis.root" );
  writer.write( vobs );

  // The End:
  return;

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

  vector<Double_t> eminFraction( 9 );
  eminFraction[0]= 0.02;
  eminFraction[1]= 0.04;
  eminFraction[2]= 0.06;
  eminFraction[3]= 0.08;
  eminFraction[4]= 0.10;
  eminFraction[5]= 0.12;
  eminFraction[6]= 0.14;
  eminFraction[7]= 0.16;
  eminFraction[8]= 0.18;
  ObsFastJetEmin* obsakteminR3= new ObsFastJetEmin( "antikteminR3", "eeantikt", 3, 0.7, 
						    eminFraction, allAnalyses );
  vector<Double_t> Rvalues( 8 );
  Rvalues[0]= 0.2;
  Rvalues[1]= 0.4;
  Rvalues[2]= 0.6;
  Rvalues[3]= 0.7;
  Rvalues[4]= 0.8;
  Rvalues[5]= 1.0;
  Rvalues[6]= 1.2;
  Rvalues[7]= 1.4;
  ObsFastJetR* obsaktRR3= new ObsFastJetR( "antiktRR3", "eeantikt", 3, 0.06, 
					   Rvalues, allAnalyses );
  vector<Observable*> vobs;
  vobs.push_back( obsakteminR3 );
  vobs.push_back( obsaktRR3 );

  vector<string> selections;
  selections.push_back( "stand" );
  selections.push_back( "costt07" );
  selections.push_back( "nch7" );
  selections.push_back( "nonrad" );
  selections.push_back( "both" );
  map<string,int> Nselected;
  for( size_t isel= 0; isel < selections.size(); isel++ ) {
    Nselected[selections[isel]]= 0;
  }

  NtupleReader* ntr= new NtupleReader( pyfilename );
  Int_t nevnt= ntr->GetNumberEntries();
  for( Int_t ievnt= 0; ievnt < TMath::Min( nevnt, maxevt ); ievnt++ ) {

    if( ntr->GetEvent( ievnt ) == 0 ) {
      cout << "event " << ievnt << " not found, skip" << endl;
      continue;
    }

    for( size_t i= 0; i < allAnalyses.size(); i++ ) {
      string cuts= allAnalyses[i].getCuts();
      string mccuts= allAnalyses[i].getMccuts();
      if( ( cuts == "none" || (ntr->LEP1Selections())[cuts] ) &&
	  ( mccuts == "none" || ntr->MCNonRad() ) ) {
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

  TCanvas* canvd= new TCanvas( "canvd", "y23 Durham plots", 900, 900 );
  canvd->Divide( 2, 2 );
  canvd->cd( 1 );
  y23dmt->Draw();
  canvd->cd( 2 );
  y23dhadsel->Draw();
  canvd->cd( 3 );
  y23dmtnr->Draw();
  canvd->cd( 4 );
  y23dmthad->Draw( "box" );

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

  Analysis a( "py", "mt", "stand" );
  Analysis b( "py", "hadron", "none", "nonrad" );
  TCanvas* canva= new TCanvas( "canva", "ee anti-kt plots", 900, 900 );
  canva->Divide( 2, 2 );
  canva->cd( 1 );
  JetratePlot( obsakteminR3, a );
  canva->cd( 2 );
  JetratePlot( obsaktRR3, a );
  canva->cd( 3 );
  JetratePlot( obsakteminR3, b );
  canva->cd( 4 );
  JetratePlot( obsaktRR3, b );

  // canv->Print( "plot.pdf" );
    
  return;

}

