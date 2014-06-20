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
class DifferentialDataStructure;

class ObsDifferential : public Observable {

public:

  ObsDifferential( const string& name,
		   const vector<Double_t>& bins );
  virtual ~ObsDifferential() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation ) = 0;
  virtual void addAnalyses( const vector<Analysis>& variations );

protected:
  
  DifferentialDataStructure* getDifferentialDataStructure( DataStructure* ) const;
  void getAndFillDifferentialDataStructure( Double_t, const string&,
					    const map<string,DataStructure*>&,
					    Double_t weight=1.0 );

  vector<Double_t> binedges;

};

#endif
