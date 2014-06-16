
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Analysis.hh"
#include "Observable.hh"
#include "ObservableFactory.hh"
#include "Unfolder.hh"
#include "OutputWriter.hh"
#include "NtupleReader.hh"
#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <string>
using std::string;
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
  vector<string> obsnames;
  obsnames.push_back( "thrust" );
  obsnames.push_back( "durhamymerge23" );
  obsnames.push_back( "jadeymerge23" );
  obsnames.push_back( "antikteminR3" );
  obsnames.push_back( "antiktRR3" );
  obsnames.push_back( "sisconeeminR3" );
  obsnames.push_back( "sisconeRR3" );
  ObservableFactory obsfac;
  vector<Observable*> vobs= obsfac.createObservables( obsnames, allAnalyses );

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
