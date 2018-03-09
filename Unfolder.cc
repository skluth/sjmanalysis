
#include "Unfolder.hh"

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

Unfolder::Unfolder( const Analysis& measured, 
		    const Analysis& measuredmc, 
		    const Analysis& hadronlevel ) :
  measuredAnalysis(measured), 
  measuredMCAnalysis(measuredmc),
  hadronlevelAnalysis(hadronlevel) {}

void Unfolder::unfold( FilledObservable* obs ) const {

  // Check that all inputs are there:
  if( not obs->containsAnalysis( hadronlevelAnalysis ) ) {
    throw runtime_error( "Unfolder::unfold: hadronlevel analysis not found: " +
			 hadronlevelAnalysis.getTag() );
  }
  if( not obs->containsAnalysis( measuredAnalysis ) ) {
    throw runtime_error( "Unfolder::unfold: measured analysis not found: " +
			 measuredAnalysis.getTag() );
  }
  if( not obs->containsAnalysis( measuredMCAnalysis ) ) {
    throw runtime_error( "Unfolder::unfold: measured MC analysis not found: " +
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
  vector<Double_t> correctionFactors= divideChecked( valuesHadronlevel,
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

  // Error matrix for normalised bin-by-bin corrected
  // distributions a la OPAL PR404:
  if( obs->makeErrorMatrices() ) {
    DifferentialDataStructure* correctedDataDds=
      dynamic_cast<DifferentialDataStructure*>( correctedData );
    if( correctedDataDds ) {
      cout << "Unfolder::unfold: error matrix for " << obs->getName() << " "
	   << correctedDataAnalysis.getTag() << endl;
      correctedData->setErrorMatrix();
      MatrixDataStructure* errorMatrix= correctedData->getErrorMatrix();
      // calculateErrorMatrix( errorMatrix,
      // 			    valuesMeasured,
      // 			    correctedValues,
      // 			    correctionFactors,
      // 			    neventsCorrected );
      calculateErrorMatrix2( errorMatrix,
      			     correctedValues,
      			     correctedErrors,
      			     neventsCorrected );
    }
    else {
      throw logic_error( "Unfolder::unfold: wrong class" );
    }
  }
  
  obs->setDataStructure( correctedData, correctedDataAnalysis );

  // The End:
  return;

}

void Unfolder::calculateErrorMatrix( MatrixDataStructure* errorMatrix,
				     const vector<Double_t> & valuesMeasured,
				     const vector<Double_t> & correctedValues,
				     const vector<Double_t> & correctionFactors,
				     Double_t neventsCorrected ) const {
  size_t nbin= valuesMeasured.size()-2;
  for( size_t ibin= 1; ibin <= nbin; ibin++ ) {
    for( size_t jbin= 1; jbin <= nbin; jbin++ ) {
      Double_t cov= 0.0;
      for( size_t kbin= 1; kbin <= nbin; kbin++ ) {
	Double_t deltaik= 0.0;
	if( ibin == kbin ) deltaik= 1.0;
	Double_t deltajk= 0.0;
	if( jbin == kbin ) deltajk= 1.0;
	cov+= pow( correctionFactors[kbin], 2 )*valuesMeasured[kbin]*
	  ( neventsCorrected*deltaik - correctedValues[ibin] )*
	  ( neventsCorrected*deltajk - correctedValues[jbin] );
      }
      cov/= pow( neventsCorrected, 4 );
      errorMatrix->setElement( ibin, jbin, cov );
    }
  }
  return;
}

// Generalised OPAL method from calculation of the
// Jacobean of the correction after normalisation:
void Unfolder::calculateErrorMatrix2( MatrixDataStructure* errorMatrix,
				      const vector<Double_t> & correctedValues,
				      const vector<Double_t> & correctedErrors,
				      Double_t neventsCorrected ) const {
  size_t nbin= correctedValues.size()-2;
  Double_t sumw= 0.0;
  for( size_t ibin= 1; ibin <= nbin; ibin++ ) sumw+= correctedValues[ibin];
  for( size_t ibin= 1; ibin <= nbin; ibin++ ) {
    for( size_t jbin= 1; jbin <= nbin; jbin++ ) {
      Double_t cov= 0.0;
      for( size_t kbin= 1; kbin <= nbin; kbin++ ) {
	Double_t deltaik= 0.0;
	if( ibin == kbin ) deltaik= 1.0;
	Double_t deltajk= 0.0;
	if( jbin == kbin ) deltajk= 1.0;
	cov+= pow( correctedErrors[kbin], 2 )*
	  ( sumw*deltaik - correctedValues[ibin] )*
	  ( sumw*deltajk - correctedValues[jbin] );
      }
      cov/= pow( neventsCorrected*sumw, 2 );
      errorMatrix->setElement( ibin, jbin, cov );
    }
  }
  return;
}

