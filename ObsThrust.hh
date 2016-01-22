#ifndef OBSTHRUST_HH
#define OBSTHRUST_HH

#include "Rtypes.h"
#include "ObsDifferential.hh"
#include "Analysis.hh"

#include <vector>
using std::vector;
#include <string>
using std::string;
#include <map>
using std::map;

class NtupleReader;
class DataStructrue;
class FilledObservable;

class ObsThrust : public ObsDifferential {

public:

  ObsThrust( const vector<Double_t>& bins, 
	     const vector<Analysis>& variations,
	     const bool lprint=true );
  ~ObsThrust() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation );
  virtual vector<FilledObservable*> getFilledObservables() const;

private:

  void addAnalyses( const vector<Analysis>& variations );

  map<string,DataStructure*> weighted1;
  map<string,DataStructure*> weighted2;

};

#endif
