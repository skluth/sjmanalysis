
#include "AnalysisProcessor.hh"

#include "SjmConfigParser.hh"
#include "Analysis.hh"
#include "Observable.hh"
#include "FilledObservable.hh"
#include "ObservableFactory.hh"
#include "Unfolder.hh"
#include "OutputWriter.hh"
#include "NtupleReader.hh"
using std::vector;
using std::string;
#include <sstream>
using std::ostringstream;
#include <iostream>
using std::cout;
using std::endl;
#include <stdexcept>
using std::logic_error;


AnalysisProcessor::AnalysisProcessor( const SjmConfigParser& sjmcp ) {
  datafilenames= sjmcp.getDataFiles();
  pyfilenames= sjmcp.getSignalMCFiles();
  hwfilenames= sjmcp.getAltSignalMCFiles();
  maxevt= sjmcp.getMaxEvent();
  obsnames= sjmcp.getObservables();
}  


void AnalysisProcessor::processAnalyses( const vector<Analysis>& analyses,
					 const vector<Observable*>& vobs,
					 const string& filename ) {
  cout << "processAnalyses: file " << filename << ", analyses:" << endl;
  for( size_t i= 0; i < analyses.size(); i++ ) {
    cout << analyses[i].getTag() << endl;
  }
  NtupleReader* ntr= new NtupleReader( filename.c_str() );
  Int_t nevnt= ntr->GetNumberEntries();
  if( nevnt > maxevt ) cout << "processAnalyses: process " 
			    << maxevt << " events" << endl;
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

void AnalysisProcessor::processUnfolding( const vector<Analysis>& measuredAnalyses, 
					  const string& unfoldsource,
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

vector<FilledObservable*> 
AnalysisProcessor::getFilled( const vector<Observable*>& vobs ) {
  vector<FilledObservable*> vfobs;
  for( size_t iobs= 0; iobs < vobs.size(); iobs++ ) {
    vector<FilledObservable*> vfobspart= vobs[iobs]->getFilledObservables();
    vfobs.insert( vfobs.end(), vfobspart.begin(), vfobspart.end() );
  }
  return vfobs;
}

void AnalysisProcessor::LEP1Analysis() {

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
  ObservableFactory obsfac;
  vector<Observable*> vobs;
  try {
    vobs= obsfac.createObservables( obsnames, allAnalyses );
  }
  catch( const std::exception& e ) {
    cout << "AnalysisProcessor::LEP1Analysis: cought exception: " << e.what() << endl;
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
    for( const string& datafilename : datafilenames ) {
      processAnalyses( measuredAnalyses, vobs, datafilename );
    }
    for( const string& pyfilename : pyfilenames ) {
      processAnalyses( pyAnalyses, vobs, pyfilename );
    }
    for( const string& hwfilename : hwfilenames ) {
      processAnalyses( hwAnalyses, vobs, hwfilename );
    }
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

