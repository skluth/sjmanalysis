
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
#include "LEPNtupleReader.hh"
#include "LEP1NtupleReader.hh"
#include "LEP2NtupleReader.hh"

// #include "HepMC2Reader.hh"
#include "HepMCRootReader.hh"

#include "VectorHelpers.hh"
#include "DataStructure.hh"
#include "fastjet/Error.hh"

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
AnalysisProcessor::processAnalyses( const vector<Analysis> & analyses,
				    const vector<Observable*> & vobs,
				    std::map<std::string,int> & cutflowCounter,
				    const string & filename ) {
  cout << "processAnalyses: file " << filename << endl;
  NtupleReader* ntr= createNtupleReader( filename );
  Int_t nevents= processAnalysesNtr( analyses, vobs, cutflowCounter, ntr );
  delete ntr;
  return nevents;
}

void countCutflow( const map<string,Bool_t> & cutStatus,
		     std::map<std::string,int> & cutflowCounter ) {
  if( cutflowCounter.empty() ) {
    std::cout << "countCutflow: initialise cut flow map" << std::endl;
    for( auto const & keyValue : cutStatus ) {
      std::cout << keyValue.first << " ";
      cutflowCounter[keyValue.first]= 0;
    }
    std::cout << std::endl;
  }
  for( auto const & keyValue : cutStatus ) {
    const string cutKey= keyValue.first;
    Bool_t cutValue= keyValue.second;
    if( cutValue ) cutflowCounter[cutKey]+= 1;
  }
}

Int_t
AnalysisProcessor::processAnalysesNtr( const vector<Analysis> & analyses,
				       const vector<Observable*> & vobs,
				       std::map<std::string,int> & cutflowCounter,
				       NtupleReader* ntr ) {
  cout << "processAnalysesNtr: analyses:" << endl;
  for( const Analysis & analysis : analyses ) cout << analysis.getTag() << endl;
  if( maxevt > 0 ) cout << "processAnalyses: process "
			<< maxevt << " events" << endl;
  string ecms= sjmConfigs.getItem<string>( "General.energy" );
  Int_t ievnt= 0;
  while( ntr->GetNextEvent( maxevt ) ) {
    const map<string,Bool_t> selections= ntr->getSelections( ecms );
    const map<string,Bool_t> cutflow= ntr->getCutflow();
    countCutflow( cutflow, cutflowCounter );
    bool MCnonrad= ntr->MCNonRad();
    for( const Analysis& analysis : analyses ) {
      string cuts= analysis.getCuts();
      string mccuts= analysis.getMccuts();
      try {
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
      catch( const std::exception& e ) {
	cout << "AnalysisProcessor::processAnalysis: filling cought exception: "
	     << e.what() << endl;
	cout << ievnt << " " << analysis.getTag() << endl;
	for( const auto & keyValue : selections ) {
	  cout << keyValue.first << " (" << keyValue.second << ")";
	}
	cout << "\ncuts: " << cuts << endl;
      }
      catch( const fastjet::Error& fe ) {
	cout << "AnalysisProcessor::processAnalysis: filling cought fastjet exception: "
	     << fe.message() << ", event " << ievnt << endl;
      }
    }
    ievnt++;
  }
  return ievnt;
}

void
AnalysisProcessor::processUnfolding( const vector<Analysis>& measuredAnalyses,
				     const string& unfoldsource,
				     const vector<FilledObservable*>& vfobs ) {
  cout << "processUnfolding: unfolding for analyses:" << endl;
  Analysis hadronlevel( unfoldsource, "hadron", "none", "nonrad" );
  cout << "Hadron level: " << hadronlevel.getTag() << endl;
  vector<string> mtxObservables=
    sjmConfigs.getItem<vector<string>>( "Observables.mtxunfold" );
  if( mtxObservables.size() > 0 ) {
    cout << "Mtx unfolding for:";
    for( const string & obsname : mtxObservables ) {
      cout << " " << obsname;
    }
    cout << endl;
  }
  for( const Analysis& measured : measuredAnalyses ) {
    Analysis measuredMC( measured );
    measuredMC.setSource( unfoldsource );
    measuredMC.setBkgStatus( "none" );
    cout << "Data, MC signal: " << measured.getTag() << ", " << measuredMC.getTag() << endl;
    BbbUnfolder bbbunfolder( measured, measuredMC, hadronlevel );
    MtxUnfolder mtxunfolder( measured, measuredMC, hadronlevel );
    for( FilledObservable* fobs : vfobs ) {
      if( std::find( mtxObservables.begin(), mtxObservables.end(),
		     fobs->getName() ) != mtxObservables.end() ) {
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

  // For all data analyses in measuredAnalyses find bkg MC analyses and subtract:
  vector<Analysis> subtractedMeasuredAnalyses;

  // All analyses
  for( const Analysis & analysis : measuredAnalyses ) {

    // Analysis tags for bkgs and subtracted analysis:
    Analysis llqqAnalysis( analysis );
    llqqAnalysis.setSource( "llqq" );
    Analysis qqqqAnalysis( analysis );
    qqqqAnalysis.setSource( "qqqq" );
    Analysis eeqqAnalysis( analysis );
    eeqqAnalysis.setSource( "eeqq" );
    Analysis subtractedAnalysis( analysis );
    subtractedAnalysis.setBkgStatus( "llqq:qqqq:eeqq" );
    subtractedMeasuredAnalyses.push_back( subtractedAnalysis );

    // All observables:
    for( FilledObservable* obs : vfobs ) {

      // Check bkgs are there:
      for( const Analysis & analysisInObs : vector<Analysis> { analysis,
	    llqqAnalysis, qqqqAnalysis, eeqqAnalysis } ) {
	if( not obs->containsAnalysis( analysisInObs ) ) {
	  throw std::runtime_error( "Bkg subtraction: analysis not found: " +
				    analysisInObs.getTag() );
	}
      }

      // Get data:
      DataStructure* data= obs->getDataStructure( analysis );
      DataStructure* llqq= obs->getDataStructure( llqqAnalysis );
      DataStructure* qqqq= obs->getDataStructure( qqqqAnalysis );
      DataStructure* eeqq= obs->getDataStructure( eeqqAnalysis );
      vector<Double_t> valuesdata= data->getValues();
      vector<Double_t> valuesllqq= llqq->getValues();
      vector<Double_t> valuesqqqq= qqqq->getValues();
      vector<Double_t> valueseeqq= eeqq->getValues();
      vector<Double_t> errorsdata= data->getErrors();
      Double_t neventsdata= data->getNEvents();
      Double_t neventsllqq= llqq->getNEvents();
      Double_t neventsqqqq= qqqq->getNEvents();
      Double_t neventseeqq= eeqq->getNEvents();

      // Subtract with scaling factor variation for default "mt stand" analysis,
      // otherwise subtract once unscaled:
      auto subtractScaled= [&]( Double_t scaleFactor,
				const Analysis & subtracted ) {
	Double_t llqqweightSf= llqqweight*scaleFactor;
	Double_t qqqqweightSf= qqqqweight*scaleFactor;
	Double_t eeqqweightSf= eeqqweight*scaleFactor;
	vector<Double_t> subtractedData= ( valuesdata -
					   valuesllqq*llqqweightSf -
					   valuesqqqq*qqqqweightSf -
					   valueseeqq*eeqqweightSf );
	vector<Double_t> errorsllqq= llqq->getErrors()*llqqweightSf;
	vector<Double_t> errorsqqqq= qqqq->getErrors()*qqqqweightSf;
	vector<Double_t> errorseeqq= eeqq->getErrors()*eeqqweightSf;
	vector<Double_t> errors= sqrt( square( errorsdata ) +
				       square( errorsllqq ) +
				       square( errorsqqqq ) +
				       square( errorseeqq ) );
	Double_t subtractedNevents= ( neventsdata -
				      neventsllqq*llqqweightSf -
				      neventsqqqq*qqqqweightSf -
				      neventseeqq*eeqqweightSf );
	DataStructure* subtractedDds= data->clone();
	subtractedDds->setValues( subtractedData );
	subtractedDds->setErrors( errors );
	subtractedDds->setNEvents( subtractedNevents );
	obs->setDataStructure( subtractedDds, subtracted );
	return;
      };
      subtractScaled( 1.0, subtractedAnalysis );
      if( analysis.getTag().find( "data mt stand" ) != std::string::npos ) {
	Analysis subtractedAnalysisHi( analysis );
	subtractedAnalysisHi.setBkgStatus( "llqq:qqqq:eeqq:hi" );
	subtractedMeasuredAnalyses.push_back( subtractedAnalysisHi );
	Analysis subtractedAnalysisLo( analysis );
	subtractedAnalysisLo.setBkgStatus( "llqq:qqqq:eeqq:lo" );
	subtractedMeasuredAnalyses.push_back( subtractedAnalysisLo );
    	subtractScaled( 1.05, subtractedAnalysisHi );
	subtractScaled( 0.95, subtractedAnalysisLo );
      }

    }
  }
  return subtractedMeasuredAnalyses;
}

Double_t AnalysisProcessor::processFiles( const string & configKey,
					  const vector<Analysis> & analyses,
					  const vector<Observable*> & vobs,
					  std::map<std::string,int> & cutFlow ) {
  Double_t eventCount= 0.0;
  for( const string & filename : sjmConfigs.getFilepath( configKey ) ) {
    eventCount+= processAnalyses( analyses, vobs, cutFlow, filename );
  }
  return eventCount;
}

void AnalysisProcessor::LEPAnalysis() {

  cout << "AnalysisProcessor::LEPAnalysis: Welcome" << endl;

  // Get analysis variations from configuration:
  string LEPAnalyses= sjmConfigs.getItem<string>( "General.analyses" );
  vector<Analysis> measuredAnalyses= fillAnalyses( LEPAnalyses+".data" );
  vector<Analysis> pyAnalyses= fillAnalyses( LEPAnalyses+".signal" );
  vector<Analysis> hwAnalyses= fillAnalyses( LEPAnalyses+".altsignal" );

  vector<Analysis> allAnalyses( measuredAnalyses );
  allAnalyses.insert( allAnalyses.end(), pyAnalyses.begin(), pyAnalyses.end() );
  allAnalyses.insert( allAnalyses.end(), hwAnalyses.begin(), hwAnalyses.end() );

  // Background analyses:
  vector<Analysis> bkgllqqAnalyses;
  vector<Analysis> bkgqqqqAnalyses;
  vector<Analysis> bkgeeqqAnalyses;
  try {
    bkgllqqAnalyses= fillAnalyses( LEPAnalyses+".bkgllqq" );
    bkgqqqqAnalyses= fillAnalyses( LEPAnalyses+".bkgqqqq" );
    bkgeeqqAnalyses= fillAnalyses( LEPAnalyses+".bkgeeqq" );
    allAnalyses.insert( allAnalyses.end(), bkgllqqAnalyses.begin(), bkgllqqAnalyses.end() );
    allAnalyses.insert( allAnalyses.end(), bkgqqqqAnalyses.begin(), bkgqqqqAnalyses.end() );
    allAnalyses.insert( allAnalyses.end(), bkgeeqqAnalyses.begin(), bkgeeqqAnalyses.end() );
  }
  catch( const std::exception & e ) {
    cout << "AnalysisProcessor::LEPAnalysis: no background analyses" << endl;
  }

  // Define observables from configuration:
  cout << "AnalysisProcessor::LEPAnalysis: create observables" << endl;
  ObservableFactory obsfac( sjmConfigs );
  vector<Observable*> vobs;
  try {
    vector<string> observables=
      sjmConfigs.getItem<vector<string>>( "Observables.observable" );
    vobs= obsfac.createObservables( observables, allAnalyses );
  }
  catch( const std::exception & e ) {
    cout << "AnalysisProcessor::LEPAnalysis: create observables cought exception: "
	 << e.what() << endl;
    return;
  }

  // Add extra analyses to observables for migration matrices where needed:
  vector<Analysis> pyMatrixExtras;
  for( const Analysis & analysis : measuredAnalyses ) {
    if( analysis.getReco() == "mt" ) {
      pyMatrixExtras.push_back( Analysis( "py", "hadron", analysis.getCuts() ) );
    }
    pyMatrixExtras.push_back( Analysis( "py",
					analysis.getReco(),
					analysis.getCuts(),
					"none",
					"hadron" ) );
  }
  pyAnalyses.insert( pyAnalyses.end(), pyMatrixExtras.begin(), pyMatrixExtras.end() );
  vector<Analysis> hwMatrixExtras;
  hwMatrixExtras.push_back( Analysis( "hw", "hadron", "stand" ) );
  hwMatrixExtras.push_back( Analysis( "hw", "mt", "stand", "none", "hadron" ) );
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

  // Fill from data and MC (PYTHIA and HERWIG) ntuples:
  std::map<string,Double_t> eventCounts;
  std::map<std::string,std::map<std::string,int>> cutflowCounters;
  cutflowCounters["Data"]= std::map<std::string,int>();
  cutflowCounters["Signal"]= std::map<std::string,int>();
  cutflowCounters["AltSignal"]= std::map<std::string,int>();
  cout << "AnalysisProcessor::LEPAnalysis: fill from ntuples" << endl;
  try {
    eventCounts["Data"]= processFiles( "Data.files", measuredAnalyses, vobs,
				       cutflowCounters["Data"] );
    eventCounts["Signal"]= processFiles( "Signal.files", pyAnalyses, vobs,
					 cutflowCounters["Signal"] );
    eventCounts["AltSignal"]= processFiles( "AltSignal.files", hwAnalyses, vobs,
					    cutflowCounters["AltSignal"] );
    // And from background if present:
    if( bkgllqqAnalyses.size() > 0 and
	bkgqqqqAnalyses.size() > 0 and
	bkgeeqqAnalyses.size() > 0 ) {
      cutflowCounters["BkgWWllqq"]= std::map<std::string,int>();
      cutflowCounters["BkgWWqqqq"]= std::map<std::string,int>();
      cutflowCounters["BkgWWeeqq"]= std::map<std::string,int>();
      eventCounts["BkgWWllqq"]= processFiles( "BkgWWllqq.files", bkgllqqAnalyses, vobs,
					      cutflowCounters["BkgWWllqq"] );
      eventCounts["BkgWWqqqq"]= processFiles( "BkgWWqqqq.files", bkgqqqqAnalyses, vobs,
					      cutflowCounters["BkgWWqqqq"] );
      eventCounts["BkgWWeeqq"]= processFiles( "BkgWWeeqq.files", bkgeeqqAnalyses, vobs,
					      cutflowCounters["BkgWWeeqq"] );
    }
  }
  catch( const std::exception& e ) {
    cout << "AnalysisProcessor::LEPAnalysis: filling cought exception: "
	 << e.what() << endl;
    return;
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
    cout << "AnalysisProcessor::LEPAnalysis: no background subtraction" << endl;
  }

  // Unfolding bin-by-bin (or Matrix):
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
    measuredPyAnalyses.push_back( Analysis( "py", "tc", "stand" ) );
    processUnfolding( measuredPyAnalyses, "py", vfobs );
    vector<Analysis> measuredHwAnalyses;
    measuredHwAnalyses.push_back( Analysis( "hw", "mt", "stand" ) );
    processUnfolding( measuredHwAnalyses, "hw", vfobs );
  }
  catch( const std::exception& e ) {
    cout << "AnalysisProcessor::LEPAnalysis: unfolding cought exception: "
	 << e.what() << endl;
    return;
  }

  // Normalise and calculate stat errors with full jacobean calculation
  for( FilledObservable* fobs : vfobs ) {
    if( sjmConfigs.getItem<bool>( "General.normalise" ) ) fobs->finalise();
    fobs->Print();
  }

  // Write root objects (TH1D or TGraphErrors, and TH2D):
  OutputWriter writer( sjmConfigs.getItem<string>( "General.outfile" ) );
  writer.write( vfobs );
  writer.writeMaps( cutflowCounters );
  writer.writeMap( eventCounts, "Eventcounts" );

  // Print event counts and cut flow
  cout << "AnalysisProcessor::LEPAnalysis: Event counts:" << endl;
  for( const auto & keyValue : eventCounts ) {
    cout << keyValue.first << ": " << keyValue.second << endl;
  }
  cout << "AnalysisProcessor::LEPAnalysis: cut flows" << endl;
  for( const auto & keyValue : cutflowCounters ) {
    const string cutflowCounterKey= keyValue.first;
    cout << cutflowCounterKey << ":" << endl;
    std::map<std::string,int> cutflowCounter= keyValue.second;
    for( const auto & cut : cutflowCounter ) {
      const string cutflowKey= cut.first;
      int cutflowValue= cut.second;
      cout << cutflowKey << ": " << cutflowValue << endl;
    }
  }

  // The End:
  return;

}

void AnalysisProcessor::MCAnalysis() {

  cout << "AnalysisProcessor::MCAnalysis: Welcome" << endl;

  // Setup analyses for MC:
  string mcname= sjmConfigs.getItem<string>( "Signal.name" );
  Analysis hadronLevel( mcname + " hadron none nonrad" );
  Analysis partonLevel( mcname + " parton none nonrad" );
  vector<Analysis> analyses { hadronLevel, partonLevel };

  // Get observables from configuration:
  cout << "AnalysisProcessor::MCAnalysis: create observables" << endl;
  ObservableFactory obsfac( sjmConfigs );
  vector<Observable*> vobs;
  try {
    vector<string> observables=
      sjmConfigs.getItem<vector<string>>( "MCObservables.observable" );
    vobs= obsfac.createObservables( observables, analyses );
  }
  catch( const std::exception & e ) {
    cout << "AnalysisProcessor::MCAnalysis: create observables cought exception: "
	 << e.what() << endl;
    return;
  }

  // Fill from hepmc2 files:
  cout << "AnalysisProcessor::MCAnalysis: fill from hepmc files" << endl;
  try {
    vector<string> filenames= sjmConfigs.getItem<vector<string>>( "Signal.files" );
    string treename= sjmConfigs.getItem<string>( "Signal.treename" );
    string branchname= sjmConfigs.getItem<string>( "Signal.branchname" );
    map<string,int> cutflow;
    for( const string & filename : filenames ) {
      // HepMC2Reader hmc2r( filename );
      HepMCRootReader hmcrr( filename, treename, branchname );
      processAnalysesNtr( analyses, vobs, cutflow, &hmcrr );
    }
  }
  catch( const std::exception& e ) {
    cout << "AnalysisProcessor::MCAnalysis: filling cought exception: "
	 << e.what() << endl;
    return;
  }

  // Get FilledObservables for further processing:
  vector<FilledObservable*> vfobs= getFilled( vobs );

  // Normalise and calculate stat errors with full jacobean calculation
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

