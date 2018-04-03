#ifndef MTXUNFOLDER_HH
#define MTXUNFOLDER_HH

#include "Unfolder.hh"

class Analysis;

class MtxUnfolder : public Unfolder {

public:

  MtxUnfolder( const Analysis & measured,
	       const Analysis & measuredmc,
	       const Analysis & hadronlevel );
  virtual ~MtxUnfolder() {}
  void unfold( FilledObservable* ) const;

};

#endif
