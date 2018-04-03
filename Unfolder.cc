
#include "Unfolder.hh"

Unfolder::Unfolder( const Analysis & measured, 
		    const Analysis & measuredmc, 
		    const Analysis & hadronlevel ) :
  measuredAnalysis(measured), 
  measuredMCAnalysis(measuredmc),
  hadronlevelAnalysis(hadronlevel) {}

