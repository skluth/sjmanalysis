#ifndef OBSEEC_HH
#define OBSEEC_HH

#include "Rtypes.h"
#include "Observable.hh"
#include "Analysis.hh"

#include <iostream>
#include <vector>
#include <string>
#include <map>
using std::vector;
using std::string;
using std::map;

class NtupleReader;
class DataStructure;
class DifferentialDataStructure;

class ObsEEC : public Observable {

public:

  ObsEEC( const string& name,
	  const vector<Double_t>& bins,
	  const vector<Analysis>& variations,
	  const bool scOpt=true,
	  const bool lprint=true );
  virtual ~ObsEEC();
  virtual vector<FilledObservable*> getFilledObservables() const;
  virtual void fill( NtupleReader* ntr, const Analysis& variation );

private:  

  virtual void addAnalysis( const Analysis& );

  vector<Double_t> binedges;
  bool selfCorrelation;
  map<string,DifferentialDataStructure*> data;

};

#endif
