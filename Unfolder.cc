
#include "Unfolder.hh"
#include "FilledObservable.hh"
#include "DataStructure.hh"
#include "VectorHelpers.hh"

#include <iostream>
#include <vector>
#include <stdexcept>


Unfolder::Unfolder( const Analysis& measured, 
		    const Analysis& measuredmc, 
		    const Analysis& hadronlevel ) :
  measuredAnalysis(measured), 
  measuredMCAnalysis(measuredmc),
  hadronlevelAnalysis(hadronlevel) {}

void Unfolder::unfold( FilledObservable* obs ) const {

  // Check that all inputs are there:
  if( not obs->containsAnalysis( hadronlevelAnalysis ) ) {
    throw std::runtime_error( "Unfolder::unfold: hadronlevel analysis not found: " +
			      hadronlevelAnalysis.getTag() );
  }
  if( not obs->containsAnalysis( measuredAnalysis ) ) {
    throw std::runtime_error( "Unfolder::unfold: measured analysis not found: " +
			      measuredAnalysis.getTag() );
  }
  if( not obs->containsAnalysis( measuredMCAnalysis ) ) {
    throw std::runtime_error( "Unfolder::unfold: measured MC analysis not found: " +
			      measuredMCAnalysis.getTag() );
  }

  // Get real and simulated data:
  DataStructure* hadronlevel= obs->getDataStructure( hadronlevelAnalysis );
  DataStructure* measured= obs->getDataStructure( measuredAnalysis );
  DataStructure* measuredMC= obs->getDataStructure( measuredMCAnalysis  );
  std::vector<Double_t> valuesHadronlevel= hadronlevel->getValues();
  std::vector<Double_t> valuesMeasuredMC= measuredMC->getValues();
  std::vector<Double_t> valuesMeasured= measured->getValues();
  Double_t neventsMeasured= measured->getNEvents();
  Double_t neventsMeasuredMC= measuredMC->getNEvents();
  Double_t neventsHadronlevel= hadronlevel->getNEvents();

  // Calculate correction:
  std::vector<Double_t> correctionFactors= divideChecked( valuesHadronlevel,
							  valuesMeasuredMC );
  std::vector<Double_t> correctedValues= multiplyVectors( valuesMeasured, 
							  correctionFactors );
  std::vector<Double_t> errorsMeasured= measured->getErrors();
  std::vector<Double_t> correctedErrors= multiplyVectors( errorsMeasured, 
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

