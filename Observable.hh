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
class MatrixDataStructure;
class FilledObservable;

class Observable {

public:

  Observable( const string& );
  Observable() {}
  virtual ~Observable();
  virtual void addAnalyses( const vector<Analysis>& ) = 0;
  virtual void fill( NtupleReader*, const Analysis& ) = 0;
  virtual void print() const;
  virtual vector<FilledObservable*> getFilledObservables() const;
  string getName() const { return name; }

protected:

  void deleteDataStructures( map<string,DataStructure*>& );
  DataStructure* getDataStructure( const string&, const map<string,DataStructure*>& ) const;
  string name;
  map<string,DataStructure*> datastructures;
  map<string,MatrixDataStructure*> matrices;

};

#endif
