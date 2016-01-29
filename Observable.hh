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
  virtual void fillAllAnalyses( NtupleReader* );
  virtual vector<FilledObservable*> getFilledObservables() const = 0;
  string getName() const { return name; }
  virtual bool containsAnalysis( const Analysis& );
  virtual void addAnalyses( const vector<Analysis>& );

  void printVectorD( const string&, const vector<Double_t>& );

protected:

  virtual void addAnalysis( const Analysis& ) = 0;

  string name;
  vector<Analysis> analyses;

};

#endif
