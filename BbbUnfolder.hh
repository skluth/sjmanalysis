#ifndef BBBUNFOLDER_HH
#define BBBUNFOLDER_HH

#include "Unfolder.hh"

class Analysis;

class BbbUnfolder : public Unfolder {

public:

  BbbUnfolder( const Analysis & measured,
	       const Analysis & measuredmc,
	       const Analysis & hadronlevel );
  virtual ~BbbUnfolder() {}
  void unfold( FilledObservable* ) const;

};

#endif
