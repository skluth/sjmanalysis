
#include "MtxUnfolder.hh"

#include "Analysis.hh"
#include "FilledObservable.hh"
#include "DifferentialDataStructure.hh"
#include "MatrixDataStructure.hh" 
#include "VectorHelpers.hh"

#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <stdexcept>
using std::runtime_error;
using std::logic_error;
#include <math.h>

MtxUnfolder::MtxUnfolder( const Analysis & measured, 
			  const Analysis & measuredmc, 
			  const Analysis & hadronlevel ) :
  Unfolder( measured, measuredmc, hadronlevel ) {}

void MtxUnfolder::unfold( FilledObservable* obs ) const {

  // Check that all inputs are there:
  if( not obs->containsAnalysis( hadronlevelAnalysis ) ) {
    throw runtime_error( "MtxUnfolder::unfold: hadronlevel analysis not found: " +
			 hadronlevelAnalysis.getTag() );
  }
  if( not obs->containsAnalysis( measuredAnalysis ) ) {
    throw runtime_error( "MtxUnfolder::unfold: measured analysis not found: " +
			 measuredAnalysis.getTag() );
  }
  if( not obs->containsAnalysis( measuredMCAnalysis ) ) {
    throw runtime_error( "MtxUnfolder::unfold: measured MC analysis not found: " +
			 measuredMCAnalysis.getTag() );
  }
  Analysis migrationAnalysis( measuredMCAnalysis );
  migrationAnalysis.setReco2( "hadron" );
  if( not obs->containsAnalysis( migrationAnalysis ) ) {
    throw runtime_error( "MtxUnfolder::unfold: MC migration analysis not found: " +
			 migrationAnalysis.getTag() );
  }
  
  // Get real and simulated data:
  DataStructure* hadronlevel= obs->getDataStructure( hadronlevelAnalysis );
  DataStructure* measured= obs->getDataStructure( measuredAnalysis );
  DataStructure* measuredMC= obs->getDataStructure( measuredMCAnalysis  );
  vector<Double_t> valuesHadronlevel= hadronlevel->getValues();
  vector<Double_t> valuesMeasuredMC= measuredMC->getValues();
  vector<Double_t> valuesMeasured= measured->getValues();
  vector<Double_t> errorsMeasured= measured->getErrors();
  Double_t neventsMeasured= measured->getNEvents();
  Double_t neventsMeasuredMC= measuredMC->getNEvents();
  Double_t neventsHadronlevel= hadronlevel->getNEvents();

  // Migration matrix and non-rad efficiency:
  MatrixDataStructure* migrationMatrix= obs->getMigrationMatrix( migrationAnalysis );
  Analysis hadronlevelRadiativeAnalysis( measuredMCAnalysis );
  hadronlevelRadiativeAnalysis.setReco( "hadron" );
  DataStructure* hadronlevelRadiative=
    obs->getDataStructure( hadronlevelRadiativeAnalysis );
  vector<Double_t> valuesHadronlevelRadiative= hadronlevelRadiative->getValues();
  vector<Double_t> hadronEfficiency= divideChecked( valuesHadronlevel,
						    valuesHadronlevelRadiative );
  migrationMatrix->normaliseColumns();
  MatrixDataStructure corrMatrix= *(applyEfficiency( *migrationMatrix, hadronEfficiency ));

  // Correction and error matrix:
  vector<Double_t> correctedValues= multiply( corrMatrix, valuesMeasured );
  MatrixDataStructure* errorMatrix= similarityVector( corrMatrix, errorsMeasured );
  Double_t neventsCorrected= neventsMeasured*neventsHadronlevel/neventsMeasuredMC;
  vector<Double_t> correctedErrors;
  for( size_t i= 0; i < errorMatrix->getNdim(); i++ ) {
    correctedErrors.push_back( sqrt( errorMatrix->getElement( i, i ) ) );
  }

  cout << "MtxUnfolder::unfold: closure test" << endl;
  vector<Double_t> correctedMCValues= multiply( corrMatrix,
						valuesMeasuredMC );
  for( size_t i= 0; i < correctedMCValues.size(); i++ ) {
    cout << correctedMCValues[i] << " " << valuesHadronlevel[i] << endl;
  }

  // Create corrected distribution:
  DataStructure* correctedData= measured->clone();
  correctedData->setValues( correctedValues );
  correctedData->setErrors( correctedErrors );
  correctedData->setNEvents( neventsCorrected );
  correctedData->setErrorMatrix( errorMatrix );
  Analysis correctedDataAnalysis( measuredAnalysis );
  correctedDataAnalysis.setUnfoldSource( hadronlevelAnalysis.getSource() );
  correctedDataAnalysis.setUnfoldMethod( "mtx" );  
  obs->setDataStructure( correctedData, correctedDataAnalysis );

  // The End:
  return;

}


