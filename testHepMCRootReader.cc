
#include "gtest/gtest.h"
#include "gmock/gmock.h"
using ::testing::Return;
using ::testing::_;

#include "HepMCRootReader.hh"
#include "TLorentzVector.h"
#include <vector>

namespace hepmcreadertests {

  class HepMCRootReaderTest : public ::testing::Test {
  public:
    //HepMCRootReaderTest() : hmcr( "/mnt/scratch/skluth/sherpablackhatl_91.root" ) {}
    HepMCRootReaderTest() : hmcr( "/mnt/scratch/skluth/Sherpa91.hepmc3.root",
				  "hepmc3_tree", "hepmc3_event" ) {}
    virtual ~HepMCRootReaderTest() {}
    HepMCRootReader hmcr;
  };

  TEST_F( HepMCRootReaderTest, testGetNumberEntries ) {
    Int_t nentries= hmcr.GetNumberEntries();
    EXPECT_EQ( nentries, 500000 );
  }
 
  TEST_F( HepMCRootReaderTest, testGetEvent ) {
    bool status= hmcr.GetEvent( 0 );
    EXPECT_EQ( status, true );
    
    hmcr.printParticlesVertices();
    
  }

  TEST_F( HepMCRootReaderTest, testGetHadron ) {
    bool status= hmcr.GetEvent( 0 );
    EXPECT_EQ( status, true );
    const std::vector<TLorentzVector> vtlv= hmcr.GetLorentzVectors( "hadron" );
    EXPECT_EQ( vtlv.size(), 32u );
  }
  TEST_F( HepMCRootReaderTest, testGetParton ) {
    // bool status= hmcr.GetEvent( 0 );
    bool status= hmcr.GetEvent( 100 );
    EXPECT_EQ( status, true );
    const std::vector<TLorentzVector> vtlv= hmcr.GetLorentzVectors( "parton" );
    EXPECT_EQ( vtlv.size(), 6u );
  }
  TEST_F( HepMCRootReaderTest, testGetISR ) {
    bool status= hmcr.GetEvent( 0 );
    EXPECT_EQ( status, true );
    const std::vector<TLorentzVector> vtlv= hmcr.GetLorentzVectors( "isr" );
    EXPECT_EQ( vtlv.size(), 2u );
  }



  

} // namespace hepmcreadertests

// main must be provided:
int main( int argc, char** argv ) {
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}

