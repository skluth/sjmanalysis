
#include "MatrixDataStructure.hh"
#include "JetrateDataStructure.hh"
#include "DifferentialDataStructure.hh"
#include "ObsPartonShower.hh"
#include "Analysis.hh"
#include "NtupleReader.hh"
#include "LEP1NtupleReader.hh"
#include "FilledObservable.hh"
#include "ThrustCalculator.hh"
#include "YnmCalculator.hh"
#include "ObsFastJetDiff.hh"
#include "ObsJetrate.hh"
#include "YcutCalculator.hh"
#include "FastJetYcutCalculator.hh"
#include "FastJetEminCalculator.hh"
#include "FastJetRCalculator.hh"
#include "FastJetPxConeEminCalculator.hh"
#include "FastJetPxConeRCalculator.hh"

#include "SjmConfigParser.hh"

#include "VectorHelpers.hh"

#include "TMath.h"

#include <iostream>
#include <sstream>
using std::stringstream;
#include <string>
using std::string;
#include <algorithm>
#include <map>
using std::map;
#include <vector>
using std::vector;

#include "gtest/gtest.h"
#include "gmock/gmock.h"
using ::testing::Return;
using ::testing::_;

namespace sjmtests {

  // MatrixDataStructure fixture:
  class MatrixDataStructureTest : public ::testing::Test {
  public:
    MatrixDataStructureTest() : 
      bins { 0., 2., 4., 6., 8., 10. }, mds( bins ),
      binedges { 0.0, 1.0, 2.0 }, mds22( binedges ) {
	mds22.setElement( 1, 1, 1.0 );
	mds22.setElement( 2, 1, 2.0 );
	mds22.setElement( 1, 2, 3.0 );
	mds22.setElement( 2, 2, 4.0 );
      }
    virtual ~MatrixDataStructureTest() {}
    vector<double> bins;
    MatrixDataStructure mds;
    vector<Double_t> binedges;
    MatrixDataStructure mds22;
  };

    vector<Double_t> binedges { 0.0, 1.0, 2.0 };
    MatrixDataStructure m( binedges );
  
  // getBinedges:
  TEST_F( MatrixDataStructureTest, testgetBinedges ) {
    vector<Double_t> mybins= mds.getBinedges();
    EXPECT_EQ( mybins, bins );
  }

  // getElement:
  TEST_F( MatrixDataStructureTest, testgetElement ) {
    size_t ndim= mds.getBinedges().size();
    EXPECT_EQ( ndim, bins.size() );
    for( size_t i= 0; i < ndim; i++ ) {
      for( size_t j= 0; j < ndim; j++ ) {
	EXPECT_EQ( mds.getElement( i, j ), 0.0 );
      }
    }
  }

  // getElement exceptions:
  TEST_F( MatrixDataStructureTest, testgetElementExceptions ) {
    EXPECT_THROW( mds.getElement( -1, 1 ), std::runtime_error );
    EXPECT_THROW( mds.getElement( 1, -1 ), std::runtime_error );
    EXPECT_THROW( mds.getElement( 11, 1 ), std::runtime_error );
    EXPECT_THROW( mds.getElement( 1, 11 ), std::runtime_error );
  }

  // fill:
  TEST_F( MatrixDataStructureTest, testfill ) {
    mds.fill( 2.5, 2.5 );
    mds.fill( 0.5, 9.5 );
    mds.fill( 0.5, 9.5 );
    EXPECT_EQ( mds.getElement( 2, 2 ), 1.0 );
    EXPECT_EQ( mds.getElement( 1, 5 ), 2.0 );
  }

  TEST_F( MatrixDataStructureTest, testNormaliseColumns ) {
    mds.fill( 2.5, 2.5 );
    mds.fill( 2.5, 4.5 );
    mds.fill( 0.5, 9.5 );
    mds.fill( 4.5, 9.5 );
    mds.normaliseColumns();
    EXPECT_EQ( mds.getElement( 2, 2 ), 0.5 );
    EXPECT_EQ( mds.getElement( 2, 3 ), 0.5 );
    EXPECT_EQ( mds.getElement( 1, 5 ), 1.0 );
    EXPECT_EQ( mds.getElement( 3, 5 ), 1.0 );
  }
  
  TEST_F( MatrixDataStructureTest, testMultiply ) {
    vector<Double_t> v { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
    size_t ndim= v.size();
    for( size_t irow= 0; irow < ndim; irow++ ) {
      for( size_t icol= 0; icol < ndim; icol++ ) {
	mds.setElement( icol, irow, ndim*irow+icol+1 );
      }
    }
    vector<Double_t> result= multiply( mds, v );
    vector<Double_t> expected { 140.0, 336.0, 532.0, 728.0, 924.0, 1120.0, 1316.0 };
    EXPECT_EQ( expected, result );
  }

  TEST_F( MatrixDataStructureTest, testSimilarityVector ) {
    vector<Double_t> v { 0.0, 1.0, 2.0, 0.0 };
    MatrixDataStructure* mvmt= similarityVector( mds22, v );
    EXPECT_EQ( 17.0, mvmt->getElement( 1, 1 ) );
    EXPECT_EQ( 35.0, mvmt->getElement( 2, 1 ) );
    EXPECT_EQ( 35.0, mvmt->getElement( 1, 2 ) );
    EXPECT_EQ( 73.0, mvmt->getElement( 2, 2 ) );
  }

  TEST_F( MatrixDataStructureTest, testApplyEfficiency ) {
    vector<Double_t> veff { 0.0, 2.0, 3.0, 0.0 };
    MatrixDataStructure veffmds22= *(applyEfficiency( mds22, veff ));
    EXPECT_EQ( 2.0, veffmds22.getElement( 1, 1 ) );
    EXPECT_EQ( 4.0, veffmds22.getElement( 2, 1 ) );
    EXPECT_EQ( 9.0, veffmds22.getElement( 1, 2 ) );
    EXPECT_EQ( 12.0, veffmds22.getElement( 2, 2 ) );    
  }

  // JetrateDataStructure for 3-jet rate:
  class JetrateDataStructureTest : public ::testing::Test {
  public:
    JetrateDataStructureTest() : 
      points{ 1., 2., 3., 4., 5., 6. }, jrds( points, 3 ) {
	vector<Double_t> njets1{ 1., 2., 3., 4., 5., 6. };
	jrds.fill( njets1 );
	jrds.fill( njets1 );
	vector<Double_t> njets2{ 1., 1., 1., 1., 1., 1. };
	jrds.fill( njets2 );
      }
    virtual ~JetrateDataStructureTest() {}
    vector<double> points;
    JetrateDataStructure jrds;
  };

  // getPoints, getValues, getErrors:
  TEST_F( JetrateDataStructureTest, testgetters ) {
    vector<Double_t> mypoints= jrds.getPoints();
    EXPECT_EQ( mypoints, points );
    vector<Double_t> values= jrds.getValues();
    EXPECT_FLOAT_EQ( values[2], 2.0 );
    vector<Double_t> errors= jrds.getErrors();
    EXPECT_FLOAT_EQ( errors[2], TMath::Sqrt( 2.0 ) );
  }

  // normalise:
  TEST_F( JetrateDataStructureTest, testnormalise ) {
    jrds.normalise();
    vector<Double_t> values= jrds.getValues();
    vector<Double_t> errors= jrds.getErrors();
    EXPECT_FLOAT_EQ( values[2], 2.0/3.0 );
    EXPECT_FLOAT_EQ( errors[2], TMath::Sqrt( 2.0/3.0*(1.0-2.0/3.0)/3.0 ) );    
  }

  // normalise exceptions:
  TEST_F( JetrateDataStructureTest, testnormaliseExceptions ) {
    jrds.normalise();
    EXPECT_THROW( jrds.normalise(), std::runtime_error );  
    JetrateDataStructure localJrds;
    EXPECT_THROW( localJrds.normalise(), std::runtime_error ); 
  }
  
  // fill exceptions:
  TEST_F( JetrateDataStructureTest, testfillExceptions ) {
    vector<Double_t> njetstooshort{ 1., 1., 1., 1., 1. };
    EXPECT_THROW( jrds.fill( njetstooshort ), std::logic_error );
    vector<Double_t> njetstoolong{ 1., 1., 1., 1., 1., 1., 1. };
    EXPECT_THROW( jrds.fill( njetstoolong ), std::logic_error );
  }

  // DifferentialDataStructure:
  class DifferentialDataStructureTest : public ::testing::Test {
  public:
    DifferentialDataStructureTest() : 
      binedges{ 0., 1., 2., 3., 4., 5., 6. }, dds( binedges ) {
	dds.fill( -1.0 );
	dds.fill( 1.5 );
	dds.fill( 1.5 );
	dds.fill( 1.5 ); 
	dds.fill( 2.5 );
	dds.fill( 2.5 );
	dds.fill( 6.5 );
     }
    virtual ~DifferentialDataStructureTest() {}
    vector<double> binedges;
    DifferentialDataStructure dds;
  };

  // getBinedges, getPoints:
  TEST_F( DifferentialDataStructureTest, testbinspoints ) {
    vector<Double_t> mybinedges= dds.getBinedges();
    EXPECT_EQ( binedges, mybinedges );
    vector<Double_t> points= dds.getPoints();
    for( size_t i= 0; i < points.size(); i++ ) {
      EXPECT_EQ( (binedges[i+1]+binedges[i])/2.0, points[i] );
    }
  }

  // getValues, getErrors:
  TEST_F( DifferentialDataStructureTest, testgetters ) {
    vector<Double_t> values= dds.getValues();
    EXPECT_EQ( 8u, values.size() );
    EXPECT_FLOAT_EQ( 1.0, values[0] );
    EXPECT_FLOAT_EQ( 3.0, values[2] );
    EXPECT_FLOAT_EQ( 2.0, values[3] );
    EXPECT_FLOAT_EQ( 1.0, values[values.size()-1] );
    vector<Double_t> errors= dds.getErrors();
    EXPECT_FLOAT_EQ( errors[2], TMath::Sqrt( 3.0 ) );
    EXPECT_FLOAT_EQ( errors[3], TMath::Sqrt( 2.0 ) );
  }

  // normalise:
  TEST_F( DifferentialDataStructureTest, testnormalise ) {
    dds.normalise();
    vector<Double_t> values= dds.getValues();
    Double_t nevent= dds.getNEvents();
    EXPECT_FLOAT_EQ( 1.0, values[0] );
    EXPECT_FLOAT_EQ( 3.0/nevent, values[2] );
    EXPECT_FLOAT_EQ( 2.0/nevent, values[3] );
    EXPECT_FLOAT_EQ( 1.0, values[values.size()-1] );
  }

  // normalise exceptions:
  TEST_F( DifferentialDataStructureTest, testnormaliseExceptions ) {
    dds.normalise();
    EXPECT_THROW( dds.normalise(), std::runtime_error );  
    DifferentialDataStructure localDds;
    EXPECT_THROW( localDds.normalise(), std::runtime_error ); 
  }
  
  // Differential observables:

  class MockNtupleReader: public LEP1NtupleReader {
  public:
    MOCK_METHOD1( getThrust, Double_t( const TString& ) );
    MOCK_METHOD3( getYmerge, Double_t( const TString&, const TString&, Int_t ) );
  };

  class ObsDiffTest : public ::testing::Test {
  public:
    ObsDiffTest() : 
      thbinedges{ 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 
	0.09, 0.12, 0.15, 0.22, 0.3, 0.5 },
      ynmbinedges{ 0.00001, 3.16227766016838e-05, 0.0001, 0.000316227766016838,
	  0.001, 0.00316227766016838, 0.01, 0.0316227766016838, 0.1,
	  0.316227766016838, 1.0 },
      analysis1( "data mt stand" ), 
      analysis2( "data mt costt07" ),
      analyses{ analysis1, analysis2 },
      obst( "thrust", thbinedges, analyses, &tcalc, false ),
      y23dcalc( "durham", 2 ),
      obsy23d( "durhamymerge23", ynmbinedges, analyses, &y23dcalc, false ),
      y23jcalc( "jade", 2 ),
      obsy23j( "jadeymerge23", ynmbinedges, analyses, &y23jcalc, false ) {}
    virtual ~ObsDiffTest() {}
    vector<double> thbinedges;
    vector<double> ynmbinedges;
    Analysis analysis1;
    Analysis analysis2;
    vector<Analysis> analyses;
    ThrustCalculator tcalc;
    ObsDifferential obst;
    YnmCalculator y23dcalc;
    ObsDifferential obsy23d;
    YnmCalculator y23jcalc;
    ObsDifferential obsy23j;
    MockNtupleReader mntr;
  };

  // containsAnalysis:
  TEST_F( ObsDiffTest, testcontainsAnalysis ) {
    EXPECT_TRUE( obst.containsAnalysis( analysis1 ) );
    EXPECT_TRUE( obst.containsAnalysis( analysis2 ) );
    Analysis analysis3( "bla blo blu" );
    EXPECT_FALSE( obst.containsAnalysis( analysis3 ) );
  }
  // getName:
  TEST_F( ObsDiffTest, testgetName ) {
    EXPECT_EQ( "thrust", obst.getName() );
    EXPECT_EQ( "durhamymerge23", obsy23d.getName() );
    EXPECT_EQ( "jadeymerge23", obsy23j.getName() );
  }

  // Function object for std::find_if below:
  class nameeq {
    string name;
  public:
    nameeq( const string & namein ) : name(namein) {}
    bool operator() ( const FilledObservable* fobs ) {
      return fobs->getName() == name;
    }
  };

  // Helper to get observable values via FilledObservable class:
  vector<Double_t> getObsValues( const string & name, 
				 const Analysis & anal,
				 const Observable & obs ) {
    vector<FilledObservable*> fobs= obs.getFilledObservables();
    nameeq namepred( name );
    vector<FilledObservable*>::iterator it= std::find_if( fobs.begin(), 
							  fobs.end(), 
							  namepred );
    if( it == fobs.end() ) {
      throw std::runtime_error( "getObsValues: not found: "+name );
    }
    return (*it)->getDataStructure( anal )->getValues();
  }

  // Fill tests via mock NtupleReader:
  TEST_F( ObsDiffTest, testfill ) {
    EXPECT_CALL( mntr, getThrust(_) ).WillOnce(Return(0.25));
    obst.fill( &mntr, analysis1 );
    vector<Double_t> tvalues= getObsValues( "thrust", analysis1, obst );
    EXPECT_EQ( 1.0, tvalues[11] );
    vector<Double_t> tw1values= getObsValues( "thrustW1", analysis1, obst );
    EXPECT_EQ( 0.25, tw1values[11] );
    vector<Double_t> tw2values= getObsValues( "thrustW2", analysis1, obst );
    EXPECT_EQ( 0.0625, tw2values[11] );
    EXPECT_CALL( mntr, getYmerge(_,_,_) ).WillOnce(Return(0.11));
    obsy23d.fill( &mntr, analysis1 );
    vector<Double_t> y23dvalues= getObsValues( "durhamymerge23", analysis1, obsy23d );
    EXPECT_EQ( 1.0, y23dvalues[9] );
    EXPECT_CALL( mntr, getYmerge(_,_,_) ).WillOnce(Return(0.11));
    obsy23j.fill( &mntr, analysis1 );
    vector<Double_t> y23jvalues= getObsValues( "jadeymerge23", analysis1, obsy23j );
    EXPECT_EQ( 1.0, y23jvalues[9] );
  }

  // Helper to fill observables from real ntuple:
  void fillFromNtuple( NtupleReader* ntr, const Analysis analysis, 
		       Observable & obs, size_t nevnt=100 ) {
    for( size_t ievnt= 0; ievnt < nevnt; ievnt++ ) {
      ntr->GetEvent( ievnt );
      const map<string,Bool_t> selections= ntr->getSelections( "91.2" );
      bool MCnonrad= ntr->MCNonRad();
      string cuts= analysis.getCuts();
      string mccuts= analysis.getMccuts();
      if( ( cuts == "none" or selections.at( cuts ) ) and
          ( mccuts == "none" or MCnonrad ) ) {
	obs.fill( ntr, analysis );
      }
    }
  }

  // ObsDifferential fill tests:
  class ObsFastJetDiffTest : public ::testing::Test {
  public:
    ObsFastJetDiffTest() :
      thbinedges{ 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 
	0.09, 0.12, 0.15, 0.22, 0.3, 0.5 },
      ynmbinedges{ 0.00001, 3.16227766016838e-05, 0.0001, 0.000316227766016838,
	  0.001, 0.00316227766016838, 0.01, 0.0316227766016838, 0.1,
	  0.316227766016838, 1.0 },
      analysismt( "py mt stand" ),
      analysistc( "py tc stand" ),
      analysist( "py tracks stand" ), 
      analysisc( "py clusters stand" ), 
      analysisp( "py parton none nonrad" ), 
      analysish( "py hadron none nonrad" ),
      analyses{ analysismt, analysistc, analysist, analysisc, analysisp, analysish },
      obst( "thrust", thbinedges, analyses, &tcalc, false ),
      obsfjdd( "durhamymergefj", "eekt", ynmbinedges, analyses, false ) {
	ntr= new LEP1NtupleReader( "mc5025_1_200.root", "h10", false );
      }
    virtual ~ObsFastJetDiffTest() { delete ntr; }
    vector<Double_t> thbinedges;
    vector<Double_t> ynmbinedges;
    Analysis analysismt;
    Analysis analysistc;
    Analysis analysist;
    Analysis analysisc;
    Analysis analysisp;
    Analysis analysish;
    vector<Analysis> analyses;
    ThrustCalculator tcalc;
    ObsDifferential obst;
    ObsFastJetDiff obsfjdd;
    NtupleReader* ntr;
  };

  // Fill tests:
  TEST_F( ObsFastJetDiffTest, testthfillmt ) {
    fillFromNtuple( ntr, analysismt, obst );
    vector<Double_t> thvalues= getObsValues( "thrust", analysismt, obst );
    vector<Double_t> thexp{ 0, 1, 14, 15, 11, 11, 11, 7, 9, 3, 5, 2, 0, 0 };
    EXPECT_EQ( thexp, thvalues );
  } 
  TEST_F( ObsFastJetDiffTest, testthfilltc ) {
    fillFromNtuple( ntr, analysistc, obst );
    vector<Double_t> thvalues= getObsValues( "thrust", analysistc, obst );
    vector<Double_t> thexp{ 0, 3, 14, 13, 10, 14, 9, 8, 6, 5, 4, 3, 0, 0 };
    EXPECT_EQ( thexp, thvalues );
  } 
  TEST_F( ObsFastJetDiffTest, testthfillt ) {
    fillFromNtuple( ntr, analysist, obst );
    vector<Double_t> thvalues= getObsValues( "thrust", analysist, obst );
    vector<Double_t> thexp{ 0, 4, 15, 12, 9, 10, 13, 6, 7, 5, 7, 1, 0, 0 };
    EXPECT_EQ( thexp, thvalues );
  } 
  TEST_F( ObsFastJetDiffTest, testthfillc ) {
    fillFromNtuple( ntr, analysisc, obst );
    vector<Double_t> thvalues= getObsValues( "thrust", analysisc, obst );
    vector<Double_t> thexp{ 0, 3, 16, 12, 8, 12, 16, 6, 3, 7, 5, 1, 0, 0 };
    EXPECT_EQ( thexp, thvalues );
  } 
  TEST_F( ObsFastJetDiffTest, testthfillp ) {
    fillFromNtuple( ntr, analysisp, obst );
    vector<Double_t> thvalues= getObsValues( "thrust", analysisp, obst );
    vector<Double_t> thexp{ 2, 14, 16, 14, 7, 9, 4, 8, 11, 3, 4, 2, 0, 0 };
    EXPECT_EQ( thexp, thvalues );
  } 
  TEST_F( ObsFastJetDiffTest, testthfillh ) {
    fillFromNtuple( ntr, analysish, obst );
    vector<Double_t> thvalues= getObsValues( "thrust", analysish, obst );
    vector<Double_t> thexp{ 0, 3, 15, 16, 13, 10, 12, 3, 9, 6, 5, 2, 0, 0 };
    EXPECT_EQ( thexp, thvalues );
  } 
  TEST_F( ObsFastJetDiffTest, testfjy23dfillmt ) {
    fillFromNtuple( ntr, analysismt, obsfjdd );
    vector<Double_t> y23dfjvalues= getObsValues( "durhamymergefj23", analysismt, 
						 obsfjdd );
    vector<Double_t> y23dfjexp{ 0, 0, 0, 0, 10, 27, 21, 17, 10, 4, 0, 0 };
    EXPECT_EQ( y23dfjexp, y23dfjvalues );
  }
  TEST_F( ObsFastJetDiffTest, testfjy23dfilltc ) {
    fillFromNtuple( ntr, analysistc, obsfjdd );
    vector<Double_t> y23dfjvalues= getObsValues( "durhamymergefj23", analysistc,
						 obsfjdd );
    vector<Double_t> y23dfjexp{ 0, 0, 0, 2, 10, 23, 24, 16, 9, 5, 0, 0 };
    EXPECT_EQ( y23dfjexp, y23dfjvalues );
  }
  TEST_F( ObsFastJetDiffTest, testjy23dfillt ) {
    fillFromNtuple( ntr, analysist, obsfjdd );
    vector<Double_t> y23dfjvalues= getObsValues( "durhamymergefj23", analysist,
						 obsfjdd );
    vector<Double_t> y23dfjexp{ 0, 0, 0, 0, 13, 18, 25, 20, 11, 2, 0, 0 };
    EXPECT_EQ( y23dfjexp, y23dfjvalues );
  } 
  TEST_F( ObsFastJetDiffTest, testjy23dfillc ) {
    fillFromNtuple( ntr, analysisc, obsfjdd );
    vector<Double_t> y23dfjvalues= getObsValues( "durhamymergefj23", analysisc,
						 obsfjdd );
    vector<Double_t> y23dfjexp{ 0, 0, 0, 3, 14, 22, 24, 13, 12, 1, 0, 0 };
    EXPECT_EQ( y23dfjexp, y23dfjvalues );
  } 
  TEST_F( ObsFastJetDiffTest, testjy23dfillp ) {
    fillFromNtuple( ntr, analysisp, obsfjdd );
    vector<Double_t> y23dfjvalues= getObsValues( "durhamymergefj23", analysisp,
						 obsfjdd );
    vector<Double_t> y23dfjexp{ 1, 1, 2, 5, 12, 26, 14, 15, 14, 4, 0, 0 };
    EXPECT_EQ( y23dfjexp, y23dfjvalues );
  } 
  TEST_F( ObsFastJetDiffTest, testjy23dfillh ) {
    fillFromNtuple( ntr, analysish, obsfjdd );
    vector<Double_t> y23dfjvalues= getObsValues( "durhamymergefj23", analysish,
						 obsfjdd );
    vector<Double_t> y23dfjexp{ 0, 0, 0, 2, 14, 30, 18, 16, 11, 3, 0, 0 };
    EXPECT_EQ( y23dfjexp, y23dfjvalues );
  } 

  // Partonshower obs:
  class ObsPartonShowerTest : public ::testing::Test {
  public:
    ObsPartonShowerTest() : 
      binedgesA14{ 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
	0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1 },
      binedgesC202{ 0.355, 0.402, 0.424, 0.446, 0.468, 0.495 },
      binedgesAS{ -0.05, 0.04, 0.1, 0.16, 0.22, 0.28, 0.34, 0.43 },
      binedgesMR{ 0.0, 0.06, 0.15, 0.38, 0.69, 1.0 },
      analysis1( "py mt stand" ), 
      analysis2( "py mt costt07" ),
      analyses{ analysis1, analysis2 },
      obsps( binedgesA14, binedgesC202, binedgesAS, binedgesMR,
	     analyses, 0.0045, 0.5, false ) {
	ntr= new LEP1NtupleReader( "mc5025_1_200.root", "h10", false );
      }
    virtual ~ObsPartonShowerTest() { delete ntr; }
    vector<double> binedgesA14;
    vector<double> binedgesC202;
    vector<double> binedgesAS;
    vector<double> binedgesMR;
    Analysis analysis1;
    Analysis analysis2;
    vector<Analysis> analyses;
    ObsPartonShower obsps;
    NtupleReader* ntr;
  };

  // fill
  TEST_F( ObsPartonShowerTest, testfill ) {
    fillFromNtuple( ntr, analysis1, obsps );
    vector<Double_t> a14values= getObsValues( "a14", analysis1, obsps );
    vector<Double_t> a14exp{ 87, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0 };
    EXPECT_EQ( a14exp, a14values );
    vector<Double_t> c202values= getObsValues( "c202", analysis1, obsps );
    vector<Double_t> c202exp{ 87, 0, 0, 0, 2, 0, 0 };
    EXPECT_EQ( c202exp, c202values );
    vector<Double_t> asvalues= getObsValues( "as", analysis1, obsps );
    vector<Double_t> asexp{ 89, 0, 0, 0, 0, 0, 0, 0, 0 };
    EXPECT_EQ( asexp, asvalues );
    vector<Double_t> mrvalues= getObsValues( "mr", analysis1, obsps );
    vector<Double_t> mrexp{ 84, 0, 0, 1, 1, 3, 0 };
    EXPECT_EQ( mrexp, mrvalues );
  }

  // Jetrate tests
  class ObsJetrTest : public ::testing::Test {
  public:
    ObsJetrTest() :
      Ycutpoints{ 0.00001, 3.16227766016838e-05, 0.0001, 0.000316227766016838,
	0.001, 0.00316227766016838, 0.01, 0.0316227766016838, 0.1,
	0.316227766016838, 1.0 },
      Eminfpoints{ 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18 },
      Rpoints{ 0.2, 0.4, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4 },
      EminPxpoints{ 2, 6, 10, 14, 18, 22, 25.5 },
      RPxpoints{ 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5 },
      analysis1( "data mt stand" ), 
      analysis2( "data mt costt07" ),
      analyses{ analysis1, analysis2 },
      yccd( "durham" ),
      yccj( "jade" ),
      fjyccd( "eekt" ),
      fjyccj( "jade" ),
      fjeminc( "eeantikt", 0.7 ),
      fjrc( "eeantikt", 0.06 ),
      fjpxceminc( 0.7 ),
      fjpxcrc( 7.0 ),
      obsycutd( "durhamycut", Ycutpoints, analyses, &yccd, false ),
      obsycutj( "jadeycut", Ycutpoints, analyses, &yccj, false ),
      obsfjycutd( "durhamycutfj", Ycutpoints, analyses, &fjyccd, false ),
      obsfjycutj( "jadeycutfj", Ycutpoints, analyses, &fjyccj, false ),
      obsfjemin( "antiktemin", Eminfpoints, analyses, &fjeminc, false ),
      obsfjr( "antiktR", Rpoints, analyses, &fjrc, false ),
      obsfjpxemin( "pxconeemin", EminPxpoints, analyses, &fjpxceminc, false ),
      obsfjpxr( "pxconer", RPxpoints, analyses, &fjpxcrc, false )
      {
	ntr= new LEP1NtupleReader( "mc5025_1_200.root", "h10", false );
      }
    virtual ~ObsJetrTest() { delete ntr; }
    vector<double> Ycutpoints;
    vector<double> Eminfpoints;
    vector<double> Rpoints;
    vector<double> EminPxpoints;
    vector<double> RPxpoints;
    Analysis analysis1;
    Analysis analysis2;
    vector<Analysis> analyses;
    YcutCalculator yccd;
    YcutCalculator yccj;
    FastJetYcutCalculator fjyccd;
    FastJetYcutCalculator fjyccj;
    FastJetEminCalculator fjeminc;
    FastJetRCalculator fjrc;
    FastJetPxConeEminCalculator fjpxceminc;
    FastJetPxConeRCalculator fjpxcrc;
    ObsJetrate obsycutd;
    ObsJetrate obsycutj;
    ObsJetrate obsfjycutd;
    ObsJetrate obsfjycutj;
    ObsJetrate obsfjemin;
    ObsJetrate obsfjr;
    ObsJetrate obsfjpxemin;
    ObsJetrate obsfjpxr;
    NtupleReader* ntr;
  };

  // Durham jetrates
  TEST_F( ObsJetrTest, testfillycutd ) {
    fillFromNtuple( ntr, analysis1, obsycutd );
    vector<Double_t> n2jets= getObsValues( "durhamycutR2", analysis1, obsycutd );
    vector<Double_t> n2jetsexp{ 0, 0, 0, 0, 10, 37, 58, 75, 85, 89, 0 };
    EXPECT_EQ( n2jetsexp, n2jets );
    vector<Double_t> n3jets= getObsValues( "durhamycutR3", analysis1, obsycutd );
    vector<Double_t> n3jetsexp{ 0, 0, 0, 8, 31, 37, 31, 14, 4, 0, 0 };
    EXPECT_EQ( n3jetsexp, n3jets );
    vector<Double_t> n4jets= getObsValues( "durhamycutR4", analysis1, obsycutd );
    vector<Double_t> n4jetsexp{ 0, 0, 2, 16, 23, 12, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n4jetsexp, n4jets );
    vector<Double_t> n5jets= getObsValues( "durhamycutR5", analysis1, obsycutd );
    vector<Double_t> n5jetsexp{ 0, 1, 7, 16, 16, 3, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n5jetsexp, n5jets );
    vector<Double_t> n6jets= getObsValues( "durhamycutR6", analysis1, obsycutd );
    vector<Double_t> n6jetsexp{ 0, 1, 8, 21, 7, 0, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n6jetsexp, n6jets );
  }
  TEST_F( ObsJetrTest, testfillfjycutd ) {
    fillFromNtuple( ntr, analysis1, obsfjycutd );
    vector<Double_t> n2jets= getObsValues( "durhamycutfjR2", analysis1, obsfjycutd );
    vector<Double_t> n2jetsexp{ 0, 0, 0, 0, 10, 37, 58, 75, 85, 89, 89 };
    EXPECT_EQ( n2jetsexp, n2jets );
    vector<Double_t> n3jets= getObsValues( "durhamycutfjR3", analysis1, obsfjycutd );
    vector<Double_t> n3jetsexp{ 0, 0, 0, 8, 31, 37, 31, 14, 4, 0, 0 };
    EXPECT_EQ( n3jetsexp, n3jets );
    vector<Double_t> n4jets= getObsValues( "durhamycutfjR4", analysis1, obsfjycutd );
    vector<Double_t> n4jetsexp{ 0, 0, 2, 16, 23, 12, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n4jetsexp, n4jets );
    vector<Double_t> n5jets= getObsValues( "durhamycutfjR5", analysis1, obsfjycutd );
    // note "6" instead of "7" in 3rd ycut point
    vector<Double_t> n5jetsexp{ 0, 1, 6, 16, 16, 3, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n5jetsexp, n5jets );
    vector<Double_t> n6jets= getObsValues( "durhamycutfjR6", analysis1, obsfjycutd );
    vector<Double_t> n6jetsexp{ 89, 88, 81, 49, 9, 0, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n6jetsexp, n6jets );
  }

  // Jade jetrates
  TEST_F( ObsJetrTest, testfillycutj ) {
    fillFromNtuple( ntr, analysis1, obsycutj );
    vector<Double_t> n2jets= getObsValues( "jadeycutR2", analysis1, obsycutj );
    vector<Double_t> n2jetsexp{ 0, 0, 0, 0, 1, 10, 26, 54, 80, 89, 0 };
    EXPECT_EQ( n2jetsexp, n2jets );
    vector<Double_t> n3jets= getObsValues( "jadeycutR3", analysis1, obsycutj );
    vector<Double_t> n3jetsexp{ 0, 1, 1, 1, 6, 29, 47, 34, 9, 0, 0 };
    EXPECT_EQ( n3jetsexp, n3jets );
    vector<Double_t> n4jets= getObsValues( "jadeycutR4", analysis1, obsycutj );
    vector<Double_t> n4jetsexp{ 0, 0, 0, 2, 18, 25, 15, 1, 0, 0, 0 };
    EXPECT_EQ( n4jetsexp, n4jets );
    vector<Double_t> n5jets= getObsValues( "jadeycutR5", analysis1, obsycutj );
    vector<Double_t> n5jetsexp{ 0, 0, 0, 6, 28, 23, 1, 0, 0, 0, 0 };
    EXPECT_EQ( n5jetsexp, n5jets );
    vector<Double_t> n6jets= getObsValues( "jadeycutR6", analysis1, obsycutj );
    vector<Double_t> n6jetsexp{ 1, 1, 4, 12, 19, 2, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n6jetsexp, n6jets );
  }
  // expectations for FastJet Jade not identical to above for njet > 2
  TEST_F( ObsJetrTest, testfillfjycutj ) {
    fillFromNtuple( ntr, analysis1, obsfjycutj );
    vector<Double_t> n2jets= getObsValues( "jadeycutfjR2", analysis1, obsfjycutj );
    vector<Double_t> n2jetsexp{ 0, 0, 0, 0, 1, 10, 26, 54, 80, 89, 34 };
    EXPECT_EQ( n2jetsexp, n2jets );
    vector<Double_t> n3jets= getObsValues( "jadeycutfjR3", analysis1, obsfjycutj );
    vector<Double_t> n3jetsexp{ 0, 0, 0, 0, 4, 29, 47, 34, 9, 0, 0 };
    EXPECT_EQ( n3jetsexp, n3jets );
    vector<Double_t> n4jets= getObsValues( "jadeycutfjR4", analysis1, obsfjycutj );
    vector<Double_t> n4jetsexp{ 0, 0, 0, 1, 17, 25, 15, 1, 0, 0, 0 };
    EXPECT_EQ( n4jetsexp, n4jets );
    vector<Double_t> n5jets= getObsValues( "jadeycutfjR5", analysis1, obsfjycutj );
    vector<Double_t> n5jetsexp{ 0, 0, 0, 5, 28, 23, 1, 0, 0, 0, 0 };
    EXPECT_EQ( n5jetsexp, n5jets );
    vector<Double_t> n6jets= getObsValues( "jadeycutfjR6", analysis1, obsfjycutj );
    vector<Double_t> n6jetsexp{ 89, 89, 89, 83, 39, 2, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n6jetsexp, n6jets );
  }

  // anti-kt Emin and R
  TEST_F( ObsJetrTest, testfillEmin ) {
    fillFromNtuple( ntr, analysis1, obsfjemin );
    vector<Double_t> n2jets= getObsValues( "antikteminR2", analysis1, obsfjemin );
    vector<Double_t> n2jetsexp{ 34, 51, 59, 63, 68, 70, 73, 74, 78 };
    EXPECT_EQ( n2jetsexp, n2jets );
    vector<Double_t> n3jets= getObsValues( "antikteminR3", analysis1, obsfjemin );
    vector<Double_t> n3jetsexp{ 36, 30, 26, 24, 21, 19, 16, 15, 11 };
    EXPECT_EQ( n3jetsexp, n3jets );
    vector<Double_t> n4jets= getObsValues( "antikteminR4", analysis1, obsfjemin );
    vector<Double_t> n4jetsexp{ 12, 7, 4, 2, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n4jetsexp, n4jets );
    vector<Double_t> n5jets= getObsValues( "antikteminR5", analysis1, obsfjemin );
    vector<Double_t> n5jetsexp{ 6, 1, 0, 0, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n5jetsexp, n5jets );
    vector<Double_t> n6jets= getObsValues( "antikteminR6", analysis1, obsfjemin );
    vector<Double_t> n6jetsexp{ 1, 0, 0, 0, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n6jetsexp, n6jets );
  }
  TEST_F( ObsJetrTest, testfillR ) {
    fillFromNtuple( ntr, analysis1, obsfjr );
    vector<Double_t> n2jets= getObsValues( "antiktRR2", analysis1, obsfjr );
    vector<Double_t> n2jetsexp{ 36, 57, 57, 59, 64, 71, 78, 83 };
    EXPECT_EQ( n2jetsexp, n2jets );
    vector<Double_t> n3jets= getObsValues( "antiktRR3", analysis1, obsfjr );
    vector<Double_t> n3jetsexp{ 33, 21, 29, 26, 24, 18, 11, 6 };
    EXPECT_EQ( n3jetsexp, n3jets );
    vector<Double_t> n4jets= getObsValues( "antiktRR4", analysis1, obsfjr );
    vector<Double_t> n4jetsexp{ 18, 10, 3, 4, 1, 0, 0, 0 };
    EXPECT_EQ( n4jetsexp, n4jets );
    vector<Double_t> n5jets= getObsValues( "antiktRR5", analysis1, obsfjr );
    vector<Double_t> n5jetsexp{ 0, 1, 0, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n5jetsexp, n5jets );
    vector<Double_t> n6jets= getObsValues( "antiktRR6", analysis1, obsfjr );
    vector<Double_t> n6jetsexp{ 2, 0, 0, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n6jetsexp, n6jets );
  }

  // PXCONE Emin and R
  TEST_F( ObsJetrTest, testfillPxEmin ) {
    fillFromNtuple( ntr, analysis1, obsfjpxemin );
    vector<Double_t> n2jets= getObsValues( "pxconeeminR2", analysis1, obsfjpxemin );
    vector<Double_t> n2jetsexp{ 51, 67, 74, 79, 85, 88, 87 };
    EXPECT_EQ( n2jetsexp, n2jets );
    vector<Double_t> n3jets= getObsValues( "pxconeeminR3", analysis1, obsfjpxemin );
    vector<Double_t> n3jetsexp{ 30, 19, 15, 10, 4, 0, 0 };
    EXPECT_EQ( n3jetsexp, n3jets );
    vector<Double_t> n4jets= getObsValues( "pxconeeminR4", analysis1, obsfjpxemin );
    vector<Double_t> n4jetsexp{ 8, 3, 0, 0, 0, 0, 0 };
    EXPECT_EQ( n4jetsexp, n4jets );
  }
  TEST_F( ObsJetrTest, testfillPxR ) {
    fillFromNtuple( ntr, analysis1, obsfjpxr );
    vector<Double_t> n2jets= getObsValues( "pxconerR2", analysis1, obsfjpxr );
    vector<Double_t> n2jetsexp{ 60, 63, 71, 73, 80, 84, 89 };
    EXPECT_EQ( n2jetsexp, n2jets );
    vector<Double_t> n3jets= getObsValues( "pxconerR3", analysis1, obsfjpxr );
    vector<Double_t> n3jetsexp{ 25, 23, 16, 16, 9, 5, 0 };
    EXPECT_EQ( n3jetsexp, n3jets );
    vector<Double_t> n4jets= getObsValues( "pxconerR4", analysis1, obsfjpxr );
    vector<Double_t> n4jetsexp{ 3, 3, 2, 0, 0, 0, 0 };
    EXPECT_EQ( n4jetsexp, n4jets );
  }

  // SjmConfigParser tests:
  class SjmConfigParserTest : public ::testing::Test {
  public:
    int argc;
    const char* argv[2];
    SjmConfigParser sjmcp;
    SjmConfigParserTest() : argc(2), argv{ "program", "poconfig.cfg" }, 
      sjmcp( argc, argv ) {} 
  };

  TEST_F( SjmConfigParserTest, testgetItemString ) {
    string cfgname= sjmcp.getItem<string>( "config" );
    EXPECT_EQ( "poconfig.cfg", cfgname );
    string energy= sjmcp.getItem<string>( "General.energy" );
    EXPECT_EQ( energy, "91.2" );
  }
  TEST_F( SjmConfigParserTest, testgetItemFloat ) {
    float datalumi= sjmcp.getItem<float>( "Data.lumi" );
    EXPECT_EQ( datalumi, 1.0 );
    float smcxsec= sjmcp.getItem<float>( "Signal.xsec" );
    EXPECT_EQ( smcxsec, 2.0 );
    float asmcxsec= sjmcp.getItem<float>( "AltSignal.xsec" );
    EXPECT_EQ( asmcxsec, 3.0 );
  }
  TEST_F( SjmConfigParserTest, testgetItemInt ) {
    int maxevt= sjmcp.getItem<int>( "General.maxevt" );
    EXPECT_EQ( 999999, maxevt );
  }
  TEST_F( SjmConfigParserTest, testgetItemBool ) {
    bool normalise= sjmcp.getItem<bool>( "General.normalise" );
    EXPECT_TRUE( normalise );
  }
  TEST_F( SjmConfigParserTest, testgetItemTokens ) {
    vector<string> obs= sjmcp.getItem<vector<string>>( "Observables.observable" );
    EXPECT_EQ( obs[0], "thrust" );
    string LEPAnalyses= sjmcp.getItem<string>( "General.analyses" );
    EXPECT_EQ( LEPAnalyses, "LEP1Analyses" );
    vector<string> da= sjmcp.getItem<vector<string>>( LEPAnalyses+".data" );
    vector<string> daexp= { "data mt stand",
			    "data mt costt07",
			    "data mt nch7",
			    "data tc stand" };
    EXPECT_EQ( da, daexp );
  }
  TEST_F( SjmConfigParserTest, testgetFilepath ) {
    vector<string> datafiles= sjmcp.getFilepath( "Data.files" );
    EXPECT_EQ( datafiles[0], "/path/to/data/da91_96_200.root" );
    EXPECT_EQ( datafiles[1], "/path/to/data/da91_97_200.root" );
    vector<string> smcfiles= sjmcp.getFilepath( "Signal.files" );
    EXPECT_EQ( smcfiles[0], "/path/to/data/mc5025_1_200.root" );
    vector<string> asmcfiles= sjmcp.getFilepath( "AltSignal.files" );
    EXPECT_EQ( asmcfiles[0], "/path/to/data/mc12406_1_200.root" );
  }
  TEST_F( SjmConfigParserTest, testgetPoints ) {
    vector<double> bins= sjmcp.getPoints( "thrust" );
    vector<double> thbins = { 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 
			      0.09, 0.12, 0.15, 0.22, 0.30, 0.50 };
    EXPECT_EQ( thbins, bins );
  }
  TEST_F( SjmConfigParserTest, testgetItemDoubles ) {
    vector<double> bins= sjmcp.getItem<vector<double>>( "Points.thrust" );
    vector<double> thbins = { 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 
			      0.09, 0.12, 0.15, 0.22, 0.30, 0.50 };
    EXPECT_EQ( thbins, bins );
  }
  TEST_F( SjmConfigParserTest, testgetMissingsection ) {
    EXPECT_THROW( sjmcp.getFilepath( "BkgWWllqq.files" ), std::runtime_error );
  }

  TEST_F( SjmConfigParserTest, testgetitemTest ) {
    string str= sjmcp.getItem<string>( "General.test" );
    EXPECT_EQ( str, "a b c" );
  }
  
  // Analysis tests:
  class AnalysisTest : public ::testing::Test {
  public:
    string analysisOptions;
    Analysis analysis;
    AnalysisTest() : analysisOptions( "a b c d e f g h" ),
		     analysis( analysisOptions ) {}
  };
  TEST_F( AnalysisTest, testget ) {
    string src= analysis.getSource() ;
    EXPECT_EQ( src, "a" );
    string reco= analysis.getReco();
    EXPECT_EQ( reco, "b" );
    string cuts= analysis.getCuts();
    EXPECT_EQ( cuts, "c" );
    string mccuts= analysis.getMccuts();
    EXPECT_EQ( mccuts, "d" );
    string reco2= analysis.getReco2();
    EXPECT_EQ( reco2, "e" );
    string bkgst= analysis.getBkgStatus();
    EXPECT_EQ( bkgst, "f" );
    string unfs= analysis.getUnfoldSource();
    EXPECT_EQ( unfs, "g" );
    string unfm= analysis.getUnfoldMethod();
    EXPECT_EQ( unfm, "h" );
  }

  class TransformTest : public ::testing::Test {
  public:
    vector<double> v1 { 4, 5, 6 };
    vector<double> v2 { 1, 2, 3 };
  };
  TEST_F( TransformTest, multiplyVector ) {
    vector<double> result= multiplyVectors( v1, v2 );
    vector<double> expectation { 4, 10, 18 };
    EXPECT_EQ( expectation, result );
  }
  TEST_F( TransformTest, subtractVector ) {
    vector<double> result= subtractVectors( v1, v2 );
    vector<double> expectation { 3, 3, 3 };
    EXPECT_EQ( expectation, result );
  }
  TEST_F( TransformTest, divideChecked ) {
    vector<double> result= divideChecked( v1, v2 );
    vector<double> expectation { 4, 5.0/2.0, 2.0 };
    EXPECT_EQ( expectation, result );
  }
  TEST_F( TransformTest, divideCheckedDividebyzero ) {
    v2[1]= 0.0;
    vector<double> result= divideChecked( v1, v2 );
    vector<double> expectation { 4, 0.0, 2 };
    EXPECT_EQ( expectation, result );
  }
  TEST_F( TransformTest, divideCheckedDividebyzeroThrow ) {
    v2[1]= 0.0;
    try {
      vector<double> result= divideChecked( v1, v2, true );
      FAIL() << "Was expecting runtime_error divide by zero" << std::endl;
    }
    catch( const std::runtime_error& e ) {
      std::string txt= e.what();
      EXPECT_EQ( txt, "divide by zero" );
    }
  }
  TEST_F( TransformTest, operatorMinus ) {
    vector<double> result= v1-v2;
    vector<double> expectation { 3, 3, 3 };
    EXPECT_EQ( expectation, result );
  }
  TEST_F( TransformTest, operatorPlus ) {
    vector<double> result= v1+v2;
    vector<double> expectation { 5, 7, 9 };
    EXPECT_EQ( expectation, result );
  }
  TEST_F( TransformTest, operatorMultiplyVectors ) {
    vector<double> result= v1*v2;
    vector<double> expectation { 4, 10, 18 };
    EXPECT_EQ( expectation, result );
  }
  TEST_F( TransformTest, operatorMultiplyVectorWeight ) {
    vector<double> expectation { 8, 10, 12 }; 
    EXPECT_EQ( expectation, v1*2.0 ); 
    EXPECT_EQ( expectation, 2.0*v1 ); 
  }
  TEST_F( TransformTest, operatorDivideVectors ) {
    vector<double> result= v1/v2;
    vector<double> expectation { 4, 5.0/2.0, 2 };
    EXPECT_EQ( expectation, result );
  }
  TEST_F( TransformTest, SqrtVector ) {
    vector<double> expectation { 2, sqrt(5.0), sqrt(6) };
    vector<double> result= sqrt( v1 );
    for( size_t i= 0; i < expectation.size(); i++ ) {
      EXPECT_FLOAT_EQ( expectation[i], result[i] );
    }
  }
  TEST_F( TransformTest, SquareVector ) {
    vector<double> expectation { 16, 25, 36 };
    vector<double> result= square( v1 );
    for( size_t i= 0; i < expectation.size(); i++ ) {
      EXPECT_FLOAT_EQ( expectation[i], result[i] );
    }
  }
  TEST_F( TransformTest, chain ) {
    vector<double> expectation { sqrt(11.0), 4.0, sqrt(21.0) };
    vector<double> result= sqrt( v1*2.0 + v2*3.0 );
    for( size_t i= 0; i < expectation.size(); i++ ) {
      EXPECT_FLOAT_EQ( expectation[i], result[i] );
    }
  }
   


  
  
} // namespace

// main must be provided:
int main( int argc, char** argv ) {
  ::testing::InitGoogleMock( &argc, argv );
  return RUN_ALL_TESTS();
}

