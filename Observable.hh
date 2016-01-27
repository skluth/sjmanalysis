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
  virtual void fill( NtupleReader*, const Analysis& ) = 0;
  //virtual void print() const;
  virtual vector<FilledObservable*> getFilledObservables() const = 0;
  string getName() const { return name; }
  virtual bool containsAnalysis( const Analysis& ) = 0;
  virtual void addAnalyses( const vector<Analysis>& ) = 0;


  void printVectorD( const string&, const vector<Double_t>& );

protected:

  DataStructure* getDataStructure( const string&, 
				   const map<string,DataStructure*>& ) const;
  void deleteDataStructures( map<string,DataStructure*>& );
  bool containsAnalysisInDataStructure( const Analysis&,
					const map<string,DataStructure*>& );

  string name;

};

#endif
