#ifndef UNFOLDER_HH
#define UNFOLDER_HH

#include "Analysis.hh"
//#include <vector>
//using std::vector;

class FilledObservable;

class Unfolder {

public:

  Unfolder( const Analysis& measured, const Analysis& measuredmc,
	    const Analysis& hadronlevel );
  ~Unfolder() {}
  void unfold( FilledObservable* ) const;

private:

  Analysis measuredAnalysis;
  Analysis measuredMCAnalysis;
  Analysis hadronlevelAnalysis;

};

#endif
