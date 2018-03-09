#ifndef OBSERVABLE_HH
#define OBSERVABLE_HH

#include "Rtypes.h"
#include "Analysis.hh"
#include <vector>
#include <string>

class NtupleReader;
// class DataStructure;
// class MatrixDataStructure;
class FilledObservable;

class Observable {

public:

  Observable( const std::string & );
  Observable() {}
  virtual ~Observable();
  virtual void fill( NtupleReader*, const Analysis & ) = 0;
  virtual void fillAllAnalyses( NtupleReader* );
  virtual std::vector<FilledObservable*> getFilledObservables() const = 0;
  std::string getName() const { return name; }
  virtual bool containsAnalysis( const Analysis & );
  virtual void addAnalyses( const std::vector<Analysis> & );

  void printVectorD( const std::string &, const std::vector<Double_t> & );

protected:

  virtual void addAnalysis( const Analysis & ) = 0;

  std::string name;
  std::vector<Analysis> analyses;

};

#endif
