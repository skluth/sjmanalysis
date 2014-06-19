#ifndef FILLEDOBSERVABLE_HH
#define FILLEDOBSERVABLE_HH

#include "Rtypes.h"
#include "Analysis.hh"
#include <string>
using std::string;
#include <map>
using std::map;

class DataStructure;
//class MatrixDataStructure;
#include "MatrixDataStructure.hh"


class FilledObservable {

public:

  FilledObservable( const string&, 
		    const map<string,DataStructure*>&,
		    const map<string,MatrixDataStructure*>& mds= map<string,MatrixDataStructure*>() );
  ~FilledObservable() {}
  void finalise();
  void print() const;
  const map<string,DataStructure*>& getData() const;
  const map<string,MatrixDataStructure*>& getMatrices() const;
  DataStructure* getDataStructure( const Analysis& ) const;
  void setDataStructure( DataStructure*, const Analysis&  );
  string getName() const { return name; }

private:

  string name;
  map<string,DataStructure*> datastructures;
  map<string,MatrixDataStructure*> matrices;

};

#endif
