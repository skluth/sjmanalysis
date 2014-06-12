#ifndef OBSERVABLE_HH
#define OBSERVABLE_HH

#include "Rtypes.h"
#include "Analysis.hh"

#include <vector>
#include <string>
#include <map>
using std::vector;
using std::string;
using std::map;

class NtupleReader;
class DataStructure;

class Observable {

public:

  Observable( string namein );
  Observable() {}
  ~Observable() {}
  virtual void addAnalyses( const vector<Analysis>& variations, 
			    const vector<Double_t>& points ) = 0;
  virtual void fill( NtupleReader* ntr, const Analysis& variation ) = 0;
  void finalise();
  void print();
  virtual map<string,DataStructure*> getData();
  virtual DataStructure* getDataStructure( const Analysis& );
  virtual void setDataStructure( DataStructure*, const Analysis&  );
  string getName();

protected:

  string name;
  map<string,DataStructure*> datastructures;

};


#endif
