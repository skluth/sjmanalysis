#ifndef OBSERVABLE_HH
#define OBSERVABLE_HH

#include "Rtypes.h"
#include "Analysis.hh"
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <map>
using std::map;

class NtupleReader;
class DataStructure;
class FilledObservable;

class Observable {

public:

  Observable( const string& );
  Observable() {}
  ~Observable() {}
  virtual void addAnalyses( const vector<Analysis>& ) = 0;
  virtual void fill( NtupleReader*, const Analysis& ) = 0;
  virtual void print() const;
  virtual vector<FilledObservable*> getFilledObservables() const;
  string getName() const { return name; }

protected:

  string name;
  map<string,DataStructure*> datastructures;

};

#endif
