#ifndef OBSDIFFERENTIAL_HH
#define OBSDIFFERENTIAL_HH

#include "Rtypes.h"
#include "Observable.hh"
#include "Analysis.hh"

#include <vector>
#include <string>
#include <map>
using std::vector;
using std::string;
using std::map;

class NtupleReader;
class DataStructure;

class ObsDifferential : public Observable {

public:

  ObsDifferential( const string& name,
		   const vector<Double_t>& bins, 
		   const vector<Analysis>& variations );
  ~ObsDifferential() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation ) = 0;
  void addAnalyses( const vector<Analysis>& variations, const vector<Double_t>& bins );

protected:
  
  void getAndFillDifferentialDataStructure( Double_t value, const string& tag );
  vector<Double_t> binedges;

};

#endif
