
#include "MatrixDataStructure.hh"
#include "JetrateDataStructure.hh"
#include "DifferentialDataStructure.hh"
#include "ObsThrust.hh"
#include "ObsPartonShower.hh"
#include "Analysis.hh"
#include "NtupleReader.hh"
#include "FilledObservable.hh"

#include "TMath.h"

#include <iostream>
#include <sstream>
using std::stringstream;
#include <string>
using std::string;
#include <algorithm>

#include "gtest/gtest.h"
#include "gmock/gmock.h"
using ::testing::Return;
using ::testing::_;

namespace sjmtests {

  // MatrixDataStructure fixture:
  class MatrixDataStructureTest : public ::testing::Test {
  public:
    MatrixDataStructureTest() : 
      bins{ 0., 2., 4., 6., 8., 10. }, mds( bins ) {}
    virtual ~MatrixDataStructureTest() {}
    vector<double> bins;
    MatrixDataStructure mds;
  };

  // getBinedges:
  TEST_F( MatrixDataStructureTest, testgetBinedges ) {
    vector<Double_t> mybins= mds.getBinedges();
    EXPECT_EQ( mybins, bins );
  }

  // getElement:
  TEST_F( MatrixDataStructureTest, testgetElement ) {
    size_t ndim= mds.getBinedges().size();
    for( size_t i= 0; i < ndim; i++ ) {
      for( size_t j= 0; j < ndim; j++ ) {
	EXPECT_EQ( mds.getElement( i, j ), 0.0 );
      }
    }
  }

  // getElement exceptions
  TEST_F( MatrixDataStructureTest, testgetElementExceptions ) {
    EXPECT_THROW( mds.getElement( -1, 1 ), std::logic_error );
    EXPECT_THROW( mds.getElement( 1, -1 ), std::logic_error );
    EXPECT_THROW( mds.getElement( 11, 1 ), std::logic_error );
    EXPECT_THROW( mds.getElement( 1, 11 ), std::logic_error );
  }

  // fill:
  TEST_F( MatrixDataStructureTest, testfill ) {
    mds.fill( 2.5, 2.5 );
    mds.fill( 0.5, 9.5 );
    mds.fill( 0.5, 9.5 );
    EXPECT_EQ( mds.getElement( 2, 2 ), 1.0 );
    EXPECT_EQ( mds.getElement( 1, 5 ), 2.0 );
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

  // getPoints, getValues, getErrors
  TEST_F( JetrateDataStructureTest, testgetters ) {
    vector<Double_t> mypoints= jrds.getPoints();
    EXPECT_EQ( points, mypoints );
    vector<Double_t> values= jrds.getValues();
    EXPECT_FLOAT_EQ( values[2], 2.0 );
    vector<Double_t> errors= jrds.getErrors();
    EXPECT_FLOAT_EQ( errors[2], TMath::Sqrt( 2.0 ) );
  }

  // normalise
  TEST_F( JetrateDataStructureTest, testnormalise ) {
    jrds.normalise();
    vector<Double_t> values= jrds.getValues();
    vector<Double_t> errors= jrds.getErrors();
    EXPECT_FLOAT_EQ( values[2], 2.0/3.0 );
    EXPECT_FLOAT_EQ( errors[2], TMath::Sqrt( 2.0/3.0*(1.0-2.0/3.0)/3.0 ) );    
  }

  // fill exceptions
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
    EXPECT_EQ( 8, values.size() );
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

  // Differential observables with thrust

  class MockNtupleReader: public NtupleReader {
  public:
    MOCK_METHOD1( getThrust, Double_t( const TString& ) );
  };

  class ObsThrustTest : public ::testing::Test {
  public:
    ObsThrustTest() : 
      binedges{ 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 
	0.09, 0.12, 0.15, 0.22, 0.3, 0.5 },
      analysis1( "data", "mt", "stand" ), 
      analysis2( "data", "mt", "costt07" ),
      analyses{ analysis1, analysis2 },
      obst( binedges, analyses, false ) {}
    virtual ~ObsThrustTest() {}
    vector<double> binedges;
    Analysis analysis1;
    Analysis analysis2;
    vector<Analysis> analyses;
    ObsThrust obst;
    MockNtupleReader mntr;
  };

  // containsAnalysis
  TEST_F( ObsThrustTest, testcontainsAnalysis ) {
    EXPECT_TRUE( obst.containsAnalysis( analysis1 ) );
    EXPECT_TRUE( obst.containsAnalysis( analysis2 ) );
    Analysis analysis3( "bla", "blo", "blu" );
    EXPECT_FALSE( obst.containsAnalysis( analysis3 ) );
  }
  // getName
  TEST_F( ObsThrustTest, testgetName ) {
    string name= "thrust";
    EXPECT_EQ( name, obst.getName() );
  }

  class nameeq {
    string name;
  public:
    nameeq( const string & namein ) : name(namein) {}
    bool operator() ( const FilledObservable* fobs ) {
      return fobs->getName() == name;
    }
  };

  vector<Double_t> getObsValues( const string & name, 
				 const Analysis & anal,
				 const Observable & obs ) {
    vector<FilledObservable*> fobs= obs.getFilledObservables();
    nameeq namepred( name );
    vector<FilledObservable*>::iterator it= std::find_if( fobs.begin(), 
							  fobs.end(), 
							  namepred );
    if( it == fobs.end() ) {
      throw std::logic_error( "getObsValues: not found: "+name );
    }
    return (*it)->getDataStructure( anal )->getValues();
  }

  // fill
  TEST_F( ObsThrustTest, testfill ) {
    EXPECT_CALL( mntr, getThrust(_) ).WillOnce(Return(0.25));
    obst.fill( &mntr, analysis1 );
    vector<Double_t> tvalues= getObsValues( "thrust", analysis1, obst );
    EXPECT_EQ( 1.0, tvalues[11] );
    vector<Double_t> tw1values= getObsValues( "thrustW1", analysis1, obst );
    EXPECT_EQ( 0.25, tw1values[11] );
    vector<Double_t> tw2values= getObsValues( "thrustW2", analysis1, obst );
    EXPECT_EQ( 0.0625, tw2values[11] );
  } 

  // Partonshower obs
  class ObsPartonShowerTest : public ::testing::Test {
  public:
    ObsPartonShowerTest() : 
      binedgesA14{ 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
	0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1 },
      binedgesC202{ 0.355, 0.402, 0.424, 0.446, 0.468, 0.495 },
      binedgesAS{ -0.05, 0.04, 0.1, 0.16, 0.22, 0.28, 0.34, 0.43 },
      binedgesMR{ 0.0, 0.06, 0.15, 0.38, 0.69, 1.0 },
      analysis1( "data", "mt", "stand" ), 
      analysis2( "data", "mt", "costt07" ),
      analyses{ analysis1, analysis2 },
      obsps( binedgesA14, binedgesC202, binedgesAS, binedgesMR,
	     analyses ) {
	ntr= new NtupleReader( "mc5025_1_200.root" );
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

  void fillFromNtuple( NtupleReader* ntr, const Analysis analysis, 
		       Observable & obs, size_t nevnt=100 ) {
    for( size_t ievnt= 0; ievnt < nevnt; ievnt++ ) {
      ntr->GetEvent( ievnt );
      map<string,Bool_t> selections= ntr->LEP1Selections();
      bool MCnonrad= ntr->MCNonRad();
      string cuts= analysis.getCuts();
      string mccuts= analysis.getMccuts();
      if( ( cuts == "none" or selections[cuts] ) and
          ( mccuts == "none" or MCnonrad ) ) {
	obs.fill( ntr, analysis );
      }
    }
  }

  // fill
  TEST_F( ObsPartonShowerTest, testfill ) {
    fillFromNtuple( ntr, analysis1, obsps );
    vector<Double_t> a14values= getObsValues( "a14", analysis1, obsps );
    vector<Double_t> a14exp{ 91, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
	0, 0, 0, 0, 0, 0, 0, 0, 0 };
    EXPECT_EQ( a14exp, a14values );
    vector<Double_t> c202values= getObsValues( "c202", analysis1, obsps );
    vector<Double_t> c202exp{ 91, 0, 0, 0, 2, 0, 0 };
    EXPECT_EQ( c202exp, c202values );
    vector<Double_t> asvalues= getObsValues( "as", analysis1, obsps );
    vector<Double_t> asexp{ 93, 0, 0, 0, 0, 0, 0, 0, 0 };
    EXPECT_EQ( asexp, asvalues );
    vector<Double_t> mrvalues= getObsValues( "mr", analysis1, obsps );
    vector<Double_t> mrexp{ 88, 0, 0, 1, 1, 3, 0 };
    EXPECT_EQ( mrexp, mrvalues );
  }


} // namespace

// main must be provided:
int main( int argc, char **argv ) {
  ::testing::InitGoogleMock( &argc, argv );
  return RUN_ALL_TESTS();
}

