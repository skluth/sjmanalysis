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
		   const vector<Double_t>& bins );
  ~ObsDifferential() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation ) = 0;
  virtual void addAnalyses( const vector<Analysis>& variations );

protected:
  
  void getAndFillDifferentialDataStructure( Double_t, const string&,
					    map<string,DataStructure*>& );
  vector<Double_t> binedges;

};

#endif
