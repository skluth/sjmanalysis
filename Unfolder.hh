#ifndef UNFOLDER_HH
#define UNFOLDER_HH

#include "Analysis.hh"

class FilledObservable;

class Unfolder {

public:

  Unfolder( const Analysis & measured,
	    const Analysis & measuredmc,
	    const Analysis & hadronlevel );
  virtual ~Unfolder() {}
  virtual void unfold( FilledObservable* ) const = 0;

protected:
  
  Analysis measuredAnalysis;
  Analysis measuredMCAnalysis;
  Analysis hadronlevelAnalysis;

};

#endif
