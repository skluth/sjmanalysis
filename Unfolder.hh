#ifndef UNFOLDER_HH
#define UNFOLDER_HH

#include "Analysis.hh"
#include <vector>
using std::vector;

class Observable;

class Unfolder {

public:

  Unfolder( const Analysis& measured, const Analysis& measuredmc,
	    const Analysis& hadronlevel );
  ~Unfolder() {}
  void unfold( Observable* ) const;

private:

  Analysis measuredAnalysis;
  Analysis measuredMCAnalysis;
  Analysis hadronlevelAnalysis;

};

#endif
