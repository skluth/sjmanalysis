
#include "Unfolder.hh"
#include "FilledObservable.hh"
#include "DataStructure.hh"
#include <stdexcept>
using std::logic_error;

Unfolder::Unfolder( const Analysis& measured, 
		    const Analysis& measuredmc, 
		    const Analysis& hadronlevel ) :
  measuredAnalysis(measured), 
  measuredMCAnalysis(measuredmc),
  hadronlevelAnalysis(hadronlevel) {}

void Unfolder::unfold( FilledObservable* obs ) const {

  // Check that all inputs are there:
  if( not obs->containsAnalysis( hadronlevelAnalysis ) ) {
    throw logic_error( "Unfolder::unfold: hadronlevel analysis not found: " +
		       hadronlevelAnalysis.getTag() );
  }
  if( not obs->containsAnalysis( measuredAnalysis ) ) {
    throw logic_error( "Unfolder::unfold: measured analysis not found: " +
		       measuredAnalysis.getTag() );
  }
  if( not obs->containsAnalysis( measuredMCAnalysis ) ) {
    throw logic_error( "Unfolder::unfold: measuredMC analysis not found: " +
		       measuredMCAnalysis.getTag() );
  }

  // Get real and simulated data:
  DataStructure* hadronlevel= obs->getDataStructure( hadronlevelAnalysis );
  DataStructure* measured= obs->getDataStructure( measuredAnalysis );
  DataStructure* measuredMC= obs->getDataStructure( measuredMCAnalysis  );
  vector<Double_t> valuesHadronlevel= hadronlevel->getValues();
  vector<Double_t> valuesMeasuredMC= measuredMC->getValues();
  vector<Double_t> valuesMeasured= measured->getValues();
  Double_t neventsMeasured= measured->getNEvents();
  Double_t neventsMeasuredMC= measuredMC->getNEvents();
  Double_t neventsHadronlevel= hadronlevel->getNEvents();

  // Calculate correction:
  vector<Double_t> correctionFactors= divideVectors( valuesHadronlevel,
						     valuesMeasuredMC );
  vector<Double_t> correctedValues= multiplyVectors( valuesMeasured, 
						     correctionFactors );
  vector<Double_t> errorsMeasured= measured->getErrors();
  vector<Double_t> correctedErrors= multiplyVectors( errorsMeasured, 
						     correctionFactors );
  Double_t neventsCorrected= neventsMeasured*neventsHadronlevel/neventsMeasuredMC;

  // Create corrected distribution:
  DataStructure* correctedData= measured->clone();
  correctedData->setValues( correctedValues );
  correctedData->setErrors( correctedErrors );
  correctedData->setNEvents( neventsCorrected );
  Analysis correctedDataAnalysis( measuredAnalysis );
  correctedDataAnalysis.setUnfoldSource( hadronlevelAnalysis.getSource() );
  correctedDataAnalysis.setUnfoldMethod( "bbb" );
  obs->setDataStructure( correctedData, correctedDataAnalysis );

  // The End:
  return;

}

