#ifndef UNFOLDER_HH
#define UNFOLDER_HH

#include "Analysis.hh"

#include "Rtypes.h"

#include <vector>

class FilledObservable;
class MatrixDataStructure;

class Unfolder {

public:

  Unfolder( const Analysis & measured,
	    const Analysis & measuredmc,
	    const Analysis & hadronlevel );
  ~Unfolder() {}
  void unfold( FilledObservable* ) const;

private:

  void calculateErrorMatrix( MatrixDataStructure* errorMatrix,
			     const std::vector<Double_t> & valuesMeasured,
			     const std::vector<Double_t> & correctedValues,
			     const std::vector<Double_t> & correctionFactors,
			     Double_t neventsCorrected ) const;

  void calculateErrorMatrix2( MatrixDataStructure* errorMatrix,
			      const std::vector<Double_t> & correctedValues,
			      const std::vector<Double_t> & correctedErrors,
			      Double_t neventsCorrected ) const;
  
  Analysis measuredAnalysis;
  Analysis measuredMCAnalysis;
  Analysis hadronlevelAnalysis;

};

#endif
