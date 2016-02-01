
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Analysis.hh"
#include "Observable.hh"
#include "FilledObservable.hh"
#include "ObservableFactory.hh"
#include "Unfolder.hh"
#include "OutputWriter.hh"
#include "NtupleReader.hh"
#include <sstream>
using std::ostringstream;
#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <stdexcept>
using std::logic_error;
#endif

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
      ostringstream txt;
      txt << "processAnalyses: event not found: " << ievnt;
      throw logic_error( txt.str() );
    }
    map<string,Bool_t> selections= ntr->LEP1Selections();
    bool MCnonrad= ntr->MCNonRad();
    for( size_t ianal= 0; ianal < analyses.size(); ianal++ ) {
      Analysis analysis= analyses[ianal];
      string cuts= analysis.getCuts();
      string mccuts= analysis.getMccuts();
      if( ( cuts == "none" or selections[cuts] ) and
	  ( mccuts == "none" or MCnonrad ) ) {
	for( size_t iobs= 0; iobs < vobs.size(); iobs++ ) {
	  Observable* obs= vobs[iobs];
	  // Not all observables have all analysis variants
	  // due to filling of transfer matrices:
	  if( obs->containsAnalysis( analysis ) ) {
	    try {
	      obs->fill( ntr, analysis );
	    }
	    catch( const std::out_of_range& oor ) {
	      cout << oor.what() << " " << obs->getName() << " " 
		   << analysis.getTag() << endl;
	    }
	  }
	}
      }
    }
  }
  delete ntr;
  return;
}

void processUnfolding( const vector<Analysis>& measuredAnalyses, 
		       string unfoldsource,
		       const vector<FilledObservable*>& vobs ) {
  cout << "processUnfolding: bin-by-bin unfolding for analyses:" << endl;
  Analysis hadronlevel( unfoldsource, "hadron", "none", "nonrad" );
  cout << "Hadron level: " << hadronlevel.getTag() << endl;
  for( size_t ianal= 0; ianal < measuredAnalyses.size(); ianal++ ) {
    Analysis measured= measuredAnalyses[ianal];
    Analysis measuredMC( measured );
    measuredMC.setSource( unfoldsource );
    cout << measured.getTag() << ", " << measuredMC.getTag() << endl;
    Unfolder unfolder( measured, measuredMC, hadronlevel );
    for( size_t iobs= 0; iobs < vobs.size(); iobs++ ) {
      unfolder.unfold( vobs[iobs] );
    }
  }
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

void LEP1Analysis( Int_t maxevt=1000, 
		   const char* datafilename="da91_96_200.root",
		   const char* pyfilename="mc5025_1_200.root", 
		   const char* hwfilename="mc12406_1_200.root" ) {

  // Load libs in root before loading this macro
  // gROOT->LoadMacro("libNtupleReaderDict.so");
  // gROOT->ProcessLine(".include /home/skluth/qcd/fastjet/fastjet-3.0.6/install/include")
  // to allow ACLIC

  // Define analysis variations:
  vector<Analysis> measuredAnalyses;
  measuredAnalyses.push_back( Analysis( "data", "mt", "stand" ) );
  measuredAnalyses.push_back( Analysis( "data", "mt", "costt07" ) );
  measuredAnalyses.push_back( Analysis( "data", "mt", "nch7" ) );
  measuredAnalyses.push_back( Analysis( "data", "tc", "stand" ) );
  vector<Analysis> pyAnalyses;
  pyAnalyses.push_back( Analysis( "py", "mt", "stand" ) );
  pyAnalyses.push_back( Analysis( "py", "mt", "costt07" ) );
  pyAnalyses.push_back( Analysis( "py", "mt", "nch7" ) );
  pyAnalyses.push_back( Analysis( "py", "tc", "stand" ) );
  pyAnalyses.push_back( Analysis( "py", "hadron", "none", "nonrad" ) );
  vector<Analysis> hwAnalyses;
  hwAnalyses.push_back( Analysis( "hw", "mt", "stand" ) );
  hwAnalyses.push_back( Analysis( "hw", "hadron", "none", "nonrad" ) );
  vector<Analysis> allAnalyses( measuredAnalyses );
  allAnalyses.insert( allAnalyses.end(), pyAnalyses.begin(), pyAnalyses.end() );
  allAnalyses.insert( allAnalyses.end(), hwAnalyses.begin(), hwAnalyses.end() );

  // Define observables:
  vector<string> obsnames;
  obsnames.push_back( "thrust" );
  obsnames.push_back( "partonshower" );
  obsnames.push_back( "durhamymerge23" );
  obsnames.push_back( "jadeymerge23" );
  obsnames.push_back( "durhamymergefj" );
  obsnames.push_back( "jadeymergefj" );
  obsnames.push_back( "durhamycutfj" );
  obsnames.push_back( "jadeycutfj" );
  obsnames.push_back( "antiktemin" );
  obsnames.push_back( "antiktR" );
  obsnames.push_back( "sisconeemin" );
  obsnames.push_back( "sisconeR" );
  obsnames.push_back( "pxconeemin" );
  obsnames.push_back( "pxconeR" );
  ObservableFactory obsfac;
  vector<Observable*> vobs;
  try {
    vobs= obsfac.createObservables( obsnames, allAnalyses );
  }
  catch( const std::exception& e ) {
    cout << "Cought exception: " << e.what() << endl;
    return;
  }

  // Add extras for migration matrices where needed:
  vector<Analysis> pyMatrixExtras;
  pyMatrixExtras.push_back( Analysis( "py", "hadron", "stand", "nonrad" ) );
  pyMatrixExtras.push_back( Analysis( "py", "mt", "stand", "nonrad", "hadron" ) );
  pyAnalyses.insert( pyAnalyses.end(), pyMatrixExtras.begin(), pyMatrixExtras.end() );
  vector<Analysis> hwMatrixExtras;
  hwMatrixExtras.push_back( Analysis( "hw", "hadron", "stand", "nonrad" ) );
  hwMatrixExtras.push_back( Analysis( "hw", "mt", "stand", "nonrad", "hadron" ) );
  hwAnalyses.insert( hwAnalyses.end(), hwMatrixExtras.begin(), hwMatrixExtras.end() );
  for( size_t iobs= 0; iobs < vobs.size(); iobs++ ) {
    Observable* obs= vobs[iobs];
    if( obs->getName() == "thrust" or
	obs->getName() == "durhamymerge23" or
	obs->getName() == "jadeymerge23" or
	obs->getName() == "partonshower" ) {
      obs->addAnalyses( pyMatrixExtras );
      obs->addAnalyses( hwMatrixExtras );
    }
  }

  // Fill from data and mc (PYTHIA and HERWIG) ntuples:
  try {
    processAnalyses( measuredAnalyses, vobs, datafilename, maxevt );
    processAnalyses( pyAnalyses, vobs, pyfilename, maxevt );
    processAnalyses( hwAnalyses, vobs, hwfilename, maxevt );
  }
  catch( const std::exception& e ) {
    cout << "Cought exception: " << e.what() << endl;
    return;
  }

  // Get FilledObservables for further processing:
  vector<FilledObservable*> vfobs= getFilled( vobs );

  // Unfolding bin-by-bin:
  try {
    // PYHTHIA based:
    processUnfolding( measuredAnalyses, "py", vfobs );
    // HERWIG based for systematic:
    vector<Analysis> measuredAnalysesHw;
    measuredAnalysesHw.push_back( Analysis( "data", "mt", "stand" ) );    
    processUnfolding( measuredAnalysesHw, "hw", vfobs );
    // MC detector level with MC as cross check for PYTHIA and HERWIG:
    vector<Analysis> measuredPyAnalyses;
    measuredPyAnalyses.push_back( Analysis( "py", "mt", "stand" ) );
    measuredPyAnalyses.push_back( Analysis( "py", "mt", "costt07" ) );
    measuredPyAnalyses.push_back( Analysis( "py", "mt", "nch7" ) );
    measuredPyAnalyses.push_back( Analysis( "py", "tc", "stand" ) );
    processUnfolding( measuredPyAnalyses, "py", vfobs );
    vector<Analysis> measuredHwAnalyses;
    measuredHwAnalyses.push_back( Analysis( "hw", "mt", "stand" ) );
    processUnfolding( measuredHwAnalyses, "hw", vfobs );
  }
  catch( const std::exception& e ) {
    cout << "Cought exception: " << e.what() << endl;
    return;
  }

  // Normalise and calculate stat errors, print
  // Normalisation only during postprocessing
  for( size_t i= 0; i < vfobs.size(); i++ ) {
    //    vfobs[i]->finalise();
    vfobs[i]->print();
  }

  // Write root objects (TH1D or TGraphErrors, and TH2D):
  OutputWriter writer( "LEP1Analysis.root" );
  writer.write( vfobs );

  // The End:
  return;

}

