
#include "AnalysisProcessor.hh"

#include "SjmConfigParser.hh"
#include "Analysis.hh"
#include "Observable.hh"
#include "FilledObservable.hh"
#include "ObservableFactory.hh"
#include "BbbUnfolder.hh"
#include "MtxUnfolder.hh"
#include "OutputWriter.hh"
#include "NtupleReader.hh"
#include "LEP1NtupleReader.hh"
#include "LEP2NtupleReader.hh"
#include "VectorHelpers.hh"

#include "DataStructure.hh"

using std::map;
using std::vector;
using std::string;
#include <sstream>
using std::ostringstream;
#include <iostream>
using std::cout;
using std::endl;
#include <stdexcept>


AnalysisProcessor::AnalysisProcessor( const SjmConfigParser& sjmcp ) :
  sjmConfigs( sjmcp ), maxevt( sjmcp.getItem<int>( "General.maxevt" ) ) {
}  

NtupleReader*
AnalysisProcessor::createNtupleReader( const string& filename ) {
  string ecms= sjmConfigs.getItem<string>( "General.energy" );
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

Int_t
AnalysisProcessor::processAnalyses( const vector<Analysis>& analyses,
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
  string ecms= sjmConfigs.getItem<string>( "General.energy" );
  Int_t ievnt= 0;
  for( ; ievnt < TMath::Min( nevnt, maxevt ); ievnt++ ) {
    if( ntr->GetEvent( ievnt ) == 0 ) {
      ostringstream txt;
      txt << "processAnalyses: event not found: " << ievnt;
      throw std::runtime_error( txt.str() );
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
  return ievnt;
}

void
AnalysisProcessor::processUnfolding( const vector<Analysis>& measuredAnalyses, 
				     const string& unfoldsource,
				     const vector<FilledObservable*>& vfobs ) {
  cout << "processUnfolding: unfolding for analyses:" << endl;
  Analysis hadronlevel( unfoldsource, "hadron", "none", "nonrad" );
  cout << "Hadron level: " << hadronlevel.getTag() << endl;
  for( const Analysis& measured : measuredAnalyses ) {
    Analysis measuredMC( measured );
    measuredMC.setSource( unfoldsource );
    measuredMC.setBkgStatus( "none" );
    cout << "Data, MC signal: " << measured.getTag() << ", " << measuredMC.getTag() << endl;
    BbbUnfolder bbbunfolder( measured, measuredMC, hadronlevel );
    MtxUnfolder mtxunfolder( measured, measuredMC, hadronlevel );
    for( FilledObservable* fobs : vfobs ) {
      if( fobs->getName() == "thrust" ) {
      	mtxunfolder.unfold( fobs );
      }
      else {
      	bbbunfolder.unfold( fobs );
      }
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


vector<Analysis>
AnalysisProcessor::subtractBackground( const vector<FilledObservable*> & vfobs,
				       const vector<Analysis> & measuredAnalyses,
				       const map<string,Double_t> & eventCounts ) {

  // Calculate weights for bkg subtraction:
  Double_t lumi= sjmConfigs.getItem<float>( "Data.lumi" );
  Double_t llqqxsec= sjmConfigs.getItem<float>( "BkgWWllqq.xsec" );
  Double_t qqqqxsec= sjmConfigs.getItem<float>( "BkgWWqqqq.xsec" );
  Double_t eeqqxsec= sjmConfigs.getItem<float>( "BkgWWeeqq.xsec" );
  Double_t llqqweight= lumi*llqqxsec/eventCounts.at( "BkgWWllqq" );
  Double_t qqqqweight= lumi*qqqqxsec/eventCounts.at( "BkgWWqqqq" );
  Double_t eeqqweight= lumi*eeqqxsec/eventCounts.at( "BkgWWeeqq" );

  // For all data analyses find bkg MC analyses and subtract:
  vector<Analysis> subtractedMeasuredAnalyses;
  for( const Analysis & analysis : measuredAnalyses ) {
    Analysis llqqAnalysis( analysis );
    llqqAnalysis.setSource( "llqq" );
    Analysis qqqqAnalysis( analysis );
    llqqAnalysis.setSource( "qqqq" );
    Analysis eeqqAnalysis( analysis );
    llqqAnalysis.setSource( "eeqq" );
    Analysis subtractedDataAnalysis( analysis );
    subtractedDataAnalysis.setBkgStatus( "llqq:qqqq:eeqq" );
    subtractedMeasuredAnalyses.push_back( subtractedDataAnalysis );
    for( FilledObservable* obs : vfobs ) {
      for( const Analysis & analysisInObs : vector<Analysis> { analysis,
	    llqqAnalysis, qqqqAnalysis, eeqqAnalysis } ) {
	if( not obs->containsAnalysis( analysisInObs ) ) {
	  throw std::runtime_error( "Bkg subtraction: analysis not found: " +
				    analysisInObs.getTag() );
	}
      }
      DataStructure* data= obs->getDataStructure( analysis );
      DataStructure* llqq= obs->getDataStructure( llqqAnalysis );
      DataStructure* qqqq= obs->getDataStructure( qqqqAnalysis );
      DataStructure* eeqq= obs->getDataStructure( eeqqAnalysis );
      vector<Double_t> valuesdata= data->getValues();
      vector<Double_t> valuesllqq= llqq->getValues();
      vector<Double_t> valuesqqqq= qqqq->getValues();
      vector<Double_t> valueseeqq= eeqq->getValues();
      vector<Double_t> subtractedData= ( valuesdata -
					 valuesllqq*llqqweight -
					 valuesqqqq*qqqqweight -
					 valueseeqq*eeqqweight );
      vector<Double_t> errorsdata= data->getErrors();
      vector<Double_t> errorsllqq= llqq->getErrors()*llqqweight;
      vector<Double_t> errorsqqqq= qqqq->getErrors()*qqqqweight;
      vector<Double_t> errorseeqq= eeqq->getErrors()*eeqqweight;
      vector<Double_t> errors= sqrt( square( errorsdata ) +
				     square( errorsllqq ) +
				     square( errorsqqqq ) +
				     square( errorseeqq ) );
      Double_t neventsdata= data->getNEvents();
      Double_t neventsllqq= llqq->getNEvents();
      Double_t neventsqqqq= qqqq->getNEvents();
      Double_t neventseeqq= eeqq->getNEvents();
      Double_t subtractedNevents= ( neventsdata -
				    neventsllqq*llqqweight -
				    neventsqqqq*qqqqweight -
				    neventseeqq*eeqqweight );
      DataStructure* subtractedDds= data->clone();
      subtractedDds->setValues( subtractedData );
      subtractedDds->setErrors( errors );
      subtractedDds->setNEvents( subtractedNevents );
      obs->setDataStructure( subtractedDds, subtractedDataAnalysis );
    }
  }
  return subtractedMeasuredAnalyses;
}


void AnalysisProcessor::LEP1Analysis() {

  cout << "AnalysisProcessor::LEP1Analysis: Welcome" << endl;

  // Get analysis variations from configuration:
  vector<Analysis> measuredAnalyses= fillAnalyses( "Analyses.data" );
  vector<Analysis> pyAnalyses= fillAnalyses( "Analyses.signal" );
  vector<Analysis> hwAnalyses= fillAnalyses( "Analyses.altsignal" );
  vector<Analysis> allAnalyses( measuredAnalyses );
  allAnalyses.insert( allAnalyses.end(), pyAnalyses.begin(), pyAnalyses.end() );
  allAnalyses.insert( allAnalyses.end(), hwAnalyses.begin(), hwAnalyses.end() );

  // Background analyses:
  vector<Analysis> bkgllqqAnalyses;
  vector<Analysis> bkgqqqqAnalyses;
  vector<Analysis> bkgeeqqAnalyses;
  try {
    bkgllqqAnalyses= fillAnalyses( "Analyses.bkgllqq" );
    bkgqqqqAnalyses= fillAnalyses( "Analyses.bkgqqqq" );
    bkgeeqqAnalyses= fillAnalyses( "Analyses.bkgeeqq" );
    allAnalyses.insert( allAnalyses.end(), bkgllqqAnalyses.begin(), bkgllqqAnalyses.end() );
    allAnalyses.insert( allAnalyses.end(), bkgqqqqAnalyses.begin(), bkgqqqqAnalyses.end() );
    allAnalyses.insert( allAnalyses.end(), bkgeeqqAnalyses.begin(), bkgeeqqAnalyses.end() );
  }
  catch( const std::exception e ) {
    cout << "AnalysisProcessor::LEP1Analysis: no background analyses" << endl;
  }
  
  // Define observables from configuration:
  cout << "AnalysisProcessor::LEP1Analysis: create observables" << endl;
  ObservableFactory obsfac( sjmConfigs );
  vector<Observable*> vobs;
  try {
    vector<string> observables= 
      sjmConfigs.getItem<vector<string>>( "Observables.observable" );
    vobs= obsfac.createObservables( observables, allAnalyses );
  }
  catch( const std::exception& e ) {
    cout << "AnalysisProcessor::LEP1Analysis: create observables cought exception: "
	 << e.what() << endl;
    return;
  }

  // Add extra analyses to observables for migration matrices where needed:
  vector<Analysis> pyMatrixExtras;
  pyMatrixExtras.push_back( Analysis( "py", "hadron", "stand" ) );
  pyMatrixExtras.push_back( Analysis( "py", "hadron", "costt07" ) );
  pyMatrixExtras.push_back( Analysis( "py", "hadron", "nch7" ) );
  pyMatrixExtras.push_back( Analysis( "py", "mt", "stand", "none", "hadron" ) );
  pyMatrixExtras.push_back( Analysis( "py", "mt", "costt07", "none", "hadron" ) );
  pyMatrixExtras.push_back( Analysis( "py", "mt", "nch7", "none", "hadron" ) );
  pyMatrixExtras.push_back( Analysis( "py", "tc", "stand", "none", "hadron" ) );
  pyAnalyses.insert( pyAnalyses.end(), pyMatrixExtras.begin(), pyMatrixExtras.end() );
  vector<Analysis> hwMatrixExtras;
  hwMatrixExtras.push_back( Analysis( "hw", "hadron", "stand" ) );
  hwMatrixExtras.push_back( Analysis( "hw", "mt", "stand", "none", "hadron" ) );
  //hwMatrixExtras.push_back( Analysis( "hw", "mt", "costt07", "none", "hadron" ) );
  //hwMatrixExtras.push_back( Analysis( "hw", "mt", "nch7", "none", "hadron" ) );
  //hwMatrixExtras.push_back( Analysis( "hw", "tc", "stand", "none", "hadron" ) );
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
  std::map<string,Double_t> eventCounts;
  cout << "AnalysisProcessor::LEP1Analysis: fill from ntuples" << endl;
  try {
    eventCounts["Data"]= 0.0;
    for( const string& datafilename : sjmConfigs.getFilepath( "Data.files" ) ) {
      eventCounts["Data"]+= processAnalyses( measuredAnalyses, vobs, datafilename );
    }
    eventCounts["Signal"]= 0.0;
    for( const string& pyfilename : sjmConfigs.getFilepath( "Signal.files" ) ) {
      eventCounts["Signal"]+= processAnalyses( pyAnalyses, vobs, pyfilename );
    }
    eventCounts["AltSignal"]= 0.0;
    for( const string& hwfilename : sjmConfigs.getFilepath( "AltSignal.files" ) ) {
      eventCounts["AltSignal"]+= processAnalyses( hwAnalyses, vobs, hwfilename );
    }
    // And from background if present:
    if( bkgllqqAnalyses.size() > 0 and
	bkgqqqqAnalyses.size() > 0 and
	bkgeeqqAnalyses.size() > 0 ) {
      eventCounts["BkgWWllqq"]= 0.0;	  
      for( const string& bkgfilename : sjmConfigs.getFilepath( "BkgWWllqq.files" ) ) {
	eventCounts["BkgWWllqq"]+= processAnalyses( bkgllqqAnalyses, vobs, bkgfilename );
      }
      eventCounts["BkgWWqqqq"]= 0.0;	  
      for( const string& bkgfilename : sjmConfigs.getFilepath( "BkgWWqqqq.files" ) ) {
	eventCounts["BkgWWqqqq"]+= processAnalyses( bkgqqqqAnalyses, vobs, bkgfilename );
      }
      eventCounts["BkgWWeeqq"]= 0.0;	  
      for( const string& bkgfilename : sjmConfigs.getFilepath( "BkgWWeeqq.files" ) ) {
	eventCounts["BkgWWeeqq"]+= processAnalyses( bkgeeqqAnalyses, vobs, bkgfilename );
      }
    }
  }
  catch( const std::exception& e ) {
    cout << "AnalysisProcessor::LEP1Analysis: filling cought exception: " << e.what() << endl;
    return;
  }
  cout << "AnalysisProcessor::LEP1Analysis: Event counts:" << endl;
  for( const auto & keyValue : eventCounts ) {
    cout << keyValue.first << ": " << keyValue.second << endl;
  }
  
  // Get FilledObservables for further processing:
  vector<FilledObservable*> vfobs= getFilled( vobs );

  // Background subtraction:
  vector<Analysis> subtractedMeasuredAnalyses( measuredAnalyses );
  if( bkgllqqAnalyses.size() > 0 and
      bkgqqqqAnalyses.size() > 0 and
      bkgeeqqAnalyses.size() > 0 ) {
    subtractedMeasuredAnalyses= subtractBackground( vfobs,
						    measuredAnalyses,
						    eventCounts );
  }
  else {
    cout << "AnalysisProcessor::LEP1Analysis: no background subtraction" << endl;
  }

  // Unfolding bin-by-bin:
  try {
    // PYTHIA based:
    processUnfolding( subtractedMeasuredAnalyses, "py", vfobs );
    // HERWIG based for systematic:
    vector<Analysis> measuredAnalysesHw;
    for( const Analysis& analysis : subtractedMeasuredAnalyses ) {
      if( analysis.getTag().find( "data mt stand" ) != std::string::npos ) {
	measuredAnalysesHw.push_back( analysis );
	break;
      }
    }
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
    cout << "AnalysisProcessor::LEP1Analysis: unfolding cought exception: "
	 << e.what() << endl;
    return;
  }

  // Normalise and calculate stat errors, print
  for( FilledObservable* fobs : vfobs ) {
    if( sjmConfigs.getItem<bool>( "General.normalise" ) ) fobs->finalise();
    fobs->Print();
  }

  // Write root objects (TH1D or TGraphErrors, and TH2D):
  OutputWriter writer( sjmConfigs.getItem<string>( "General.outfile" ) );
  writer.write( vfobs );

  // The End:
  return;

}

