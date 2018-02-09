
#include "AnalysisProcessor.hh"

#include "SjmConfigParser.hh"
#include "Analysis.hh"
#include "Observable.hh"
#include "FilledObservable.hh"
#include "ObservableFactory.hh"
#include "Unfolder.hh"
#include "OutputWriter.hh"
#include "NtupleReader.hh"
#include "LEP1NtupleReader.hh"
#include "LEP2NtupleReader.hh"
using std::vector;
using std::string;
#include <sstream>
using std::ostringstream;
#include <iostream>
using std::cout;
using std::endl;
#include <stdexcept>
using std::logic_error;


AnalysisProcessor::AnalysisProcessor( const SjmConfigParser& sjmcp ) :
  sjmConfigs( sjmcp ), maxevt( sjmcp.getItem<int>( "maxevt" ) ) {
}  

NtupleReader* AnalysisProcessor::createNtupleReader( const string& filename ) {
  string ecms= sjmConfigs.getItem<string>( "energy" );
  vector<string> lep2ecms= { "130", "136", "161", "172", "183", "189", 
			     "192", "196", "200", "202", "205", "207" };
  NtupleReader* result= 0;
  if( ecms == "91.2" ) {
    cout << "AnalysisProcessor::createNtupleReader: LEP1NtupleReader" << endl;
    result= new LEP1NtupleReader( filename.c_str() );
  }
  else if( std::find( lep2ecms.begin(), lep2ecms.end(), ecms ) != lep2ecms.end() ) {
    cout << "AnalysisProcessor::createNtupleReader: LEP2NtupleReader" << endl;
    result= new LEP2NtupleReader( filename.c_str() );
  }
  else {
    throw std::runtime_error( "AnalysisProcessor::createNtupleReader: wrong ecms "+ecms );
  }
  return result;
}

void AnalysisProcessor::processAnalyses( const vector<Analysis>& analyses,
					 const vector<Observable*>& vobs,
					 const string& filename ) {
  cout << "processAnalyses: file " << filename << ", analyses:" << endl;
  for( const Analysis& analysis : analyses ) {
    cout << analysis.getTag() << endl;
  }
  NtupleReader* ntr= createNtupleReader( filename );
  Int_t nevnt= ntr->GetNumberEntries();
  if( nevnt > maxevt ) cout << "processAnalyses: process " 
			    << maxevt << " events" << endl;
  string ecms= sjmConfigs.getItem<string>( "energy" );
  for( Int_t ievnt= 0; ievnt < TMath::Min( nevnt, maxevt ); ievnt++ ) {
    if( ntr->GetEvent( ievnt ) == 0 ) {
      ostringstream txt;
      txt << "processAnalyses: event not found: " << ievnt;
      throw logic_error( txt.str() );
    }
    const map<string,Bool_t> selections= ntr->getSelections( ecms );
    bool MCnonrad= ntr->MCNonRad();
    for( const Analysis& analysis : analyses ) {
      string cuts= analysis.getCuts();
      string mccuts= analysis.getMccuts();
      if( ( cuts == "none" or selections.at( cuts ) ) and
	  ( mccuts == "none" or MCnonrad ) ) {
	for( Observable* obs : vobs ) {
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
					  const vector<FilledObservable*>& vfobs ) {
  cout << "processUnfolding: bin-by-bin unfolding for analyses:" << endl;
  Analysis hadronlevel( unfoldsource, "hadron", "none", "nonrad" );
  cout << "Hadron level: " << hadronlevel.getTag() << endl;
  for( const Analysis& measured : measuredAnalyses ) {
    Analysis measuredMC( measured );
    measuredMC.setSource( unfoldsource );
    cout << measured.getTag() << ", " << measuredMC.getTag() << endl;
    Unfolder unfolder( measured, measuredMC, hadronlevel );
    for( FilledObservable* fobs : vfobs ) {
      unfolder.unfold( fobs );
    }
  }
  return;
}

vector<FilledObservable*> 
AnalysisProcessor::getFilled( const vector<Observable*>& vobs ) {
  vector<FilledObservable*> vfobs;
  for( Observable* obs : vobs ) {
    vector<FilledObservable*> vfobspart= obs->getFilledObservables();
    vfobs.insert( vfobs.end(), vfobspart.begin(), vfobspart.end() );
  }
  return vfobs;
}

vector<Analysis> AnalysisProcessor::fillAnalyses( const string& tag ) {
  vector<Analysis> result;
  vector<string> AnalysisOptions= sjmConfigs.getItem<vector<string>>( tag );
  for( const string& options : AnalysisOptions ) {
    result.push_back( Analysis( options ) );
  }  
  return result;
}

void AnalysisProcessor::LEP1Analysis() {

  // Get analysis variations from configuration:
  vector<Analysis> measuredAnalyses= fillAnalyses( "Analyses.data" );
  vector<Analysis> pyAnalyses= fillAnalyses( "Analyses.signal" );
  vector<Analysis> hwAnalyses= fillAnalyses( "Analyses.altsignal" );
  vector<Analysis> allAnalyses( measuredAnalyses );
  allAnalyses.insert( allAnalyses.end(), pyAnalyses.begin(), pyAnalyses.end() );
  allAnalyses.insert( allAnalyses.end(), hwAnalyses.begin(), hwAnalyses.end() );

  // Define observables from configuration:
  ObservableFactory obsfac( sjmConfigs );
  vector<Observable*> vobs;
  try {
    vector<string> observables= 
      sjmConfigs.getItem<vector<string>>( "Observables.observable" );
    vobs= obsfac.createObservables( observables, allAnalyses );
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
  for( Observable* obs : vobs ) {
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
    for( const string& datafilename : sjmConfigs.getFilepath( "Data.files" ) ) {
      processAnalyses( measuredAnalyses, vobs, datafilename );
    }
    for( const string& pyfilename : sjmConfigs.getFilepath( "SignalMC.files" ) ) {
      processAnalyses( pyAnalyses, vobs, pyfilename );
    }
    for( const string& hwfilename : sjmConfigs.getFilepath( "AltSignalMC.files" ) ) {
      processAnalyses( hwAnalyses, vobs, hwfilename );
    }
  }
  catch( const std::exception& e ) {
    cout << "AnalysisProcessor::LEP1Analysis: cought exception " << e.what() << endl;
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
    // measuredPyAnalyses.push_back( Analysis( "py", "mt", "nch7" ) );
    measuredPyAnalyses.push_back( Analysis( "py", "tc", "stand" ) );
    processUnfolding( measuredPyAnalyses, "py", vfobs );
    vector<Analysis> measuredHwAnalyses;
    measuredHwAnalyses.push_back( Analysis( "hw", "mt", "stand" ) );
    processUnfolding( measuredHwAnalyses, "hw", vfobs );
  }
  catch( const std::exception& e ) {
    cout << "AnalysisProcessor::LEP1Analysis: cought exception: " << e.what() << endl;
    return;
  }

  // Normalise and calculate stat errors, print
  // Normalisation only during postprocessing
  for( FilledObservable* fobs : vfobs ) {
    //    fobs->finalise();
    fobs->print();
  }

  // Write root objects (TH1D or TGraphErrors, and TH2D):
  OutputWriter writer( sjmConfigs.getItem<string>( "outfile" ) );
  writer.write( vfobs );

  // The End:
  return;

}

