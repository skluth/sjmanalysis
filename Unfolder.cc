
#include "Unfolder.hh"
#include "FilledObservable.hh"
#include "DataStructure.hh"

Unfolder::Unfolder( const Analysis& measured, 
		    const Analysis& measuredmc, 
		    const Analysis& hadronlevel ) :
  measuredAnalysis(measured), 
  measuredMCAnalysis(measuredmc),
  hadronlevelAnalysis(hadronlevel) {}

void Unfolder::unfold( FilledObservable* obs ) const {
  DataStructure* hadronlevel= obs->getDataStructure( hadronlevelAnalysis );
  vector<Double_t> valuesHadronlevel= hadronlevel->getValues();
  DataStructure* measuredMC= obs->getDataStructure( measuredMCAnalysis  );
  vector<Double_t> valuesMeasuredMC= measuredMC->getValues();
  vector<Double_t> correctionFactors= divideVectors( valuesHadronlevel,
						     valuesMeasuredMC );
  DataStructure* measured= obs->getDataStructure( measuredAnalysis );
  vector<Double_t> valuesMeasured= measured->getValues();
  vector<Double_t> correctedValues= multiplyVectors( valuesMeasured, 
						     correctionFactors );
  DataStructure* correctedData= measured->clone();
  correctedData->setValues( correctedValues );
  vector<Double_t> errorsMeasured= measured->getErrors();
  vector<Double_t> correctedErrors= multiplyVectors( errorsMeasured, 
						     correctionFactors );
  correctedData->setErrors( correctedErrors );
  Double_t neventsMeasured= measured->getNEvents();
  Double_t neventsMeasuredMC= measuredMC->getNEvents();
  Double_t neventsHadronlevel= hadronlevel->getNEvents();
  Double_t neventsCorrected= neventsMeasured*neventsHadronlevel/neventsMeasuredMC;
  correctedData->setNEvents( neventsCorrected );
  Analysis correctedDataAnalysis( measuredAnalysis.getSource(),
				  measuredAnalysis.getReco(),
				  measuredAnalysis.getCuts(),
				  measuredAnalysis.getMccuts(),
				  hadronlevelAnalysis.getSource(),
				  "bbb" );
  obs->setDataStructure( correctedData, correctedDataAnalysis );
  return;
}

