#ifndef OBSDIFFERENTIAL_HH
#define OBSDIFFERENTIAL_HH

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
class DifferentialCalculator;

class ObsDifferential : public Observable {

public:

  ObsDifferential( const string& name,
		   const vector<Double_t>& bins,
		   const vector<Analysis>& variations,
		   const DifferentialCalculator* calc,
		   const bool lprint=true );
  virtual ~ObsDifferential();
  virtual vector<FilledObservable*> getFilledObservables() const;
  virtual bool containsAnalysis( const Analysis& );
  virtual void fill( NtupleReader* ntr, const Analysis& variation );
  virtual void addAnalyses( const vector<Analysis>& );

  //protected:
private:  

  vector<Double_t> binedges;
  map<string,DifferentialDataStructure*> data;
  map<string,DifferentialDataStructure*> weighted1;
  map<string,DifferentialDataStructure*> weighted2;
  map<string,MatrixDataStructure*> matrices;
  const DifferentialCalculator* calculator;

};

#endif
