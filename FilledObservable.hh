#ifndef FILLEDOBSERVABLE_HH
#define FILLEDOBSERVABLE_HH

#include "Rtypes.h"
#include "Analysis.hh"
#include <string>
using std::string;
#include <map>
using std::map;

class DataStructure;

class FilledObservable {

public:

  FilledObservable( const string&, const map<string,DataStructure*>& );
  FilledObservable() {}
  ~FilledObservable() {}
  void finalise();
  void print();
  map<string,DataStructure*> getData();
  DataStructure* getDataStructure( const Analysis& );
  void setDataStructure( DataStructure*, const Analysis&  );
  string getName() { return name; }

private:

  string name;
  map<string,DataStructure*> datastructures;

};

#endif
