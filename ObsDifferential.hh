#ifndef OBSDIFFERENTIAL_HH
#define OBSDIFFERENTIAL_HH

#include "Rtypes.h"
#include "Observable.hh"
#include "Analysis.hh"

#include <iostream>
#include <vector>
#include <string>
#include <map>

class NtupleReader;
class FilledObservable;
class DataStructure;
class DifferentialDataStructure;
class MatrixDataStructure;
class DifferentialCalculator;

class ObsDifferential : public Observable {

public:

  ObsDifferential( const std::string& name,
		   const std::vector<Double_t>& bins,
		   const std::vector<Analysis>& variations,
		   const DifferentialCalculator* calc,
		   const bool lprint=true );
  virtual ~ObsDifferential();
  virtual std::vector<FilledObservable*> getFilledObservables() const;
  virtual void fill( NtupleReader* ntr, const Analysis& variation );

private:  

  virtual void addAnalysis( const Analysis& );

  std::vector<Double_t> binedges;
  std::map<std::string,DifferentialDataStructure*> data;
  std::map<std::string,DifferentialDataStructure*> weighted1;
  std::map<std::string,DifferentialDataStructure*> weighted2;
  std::map<std::string,MatrixDataStructure*> matrices;
  const DifferentialCalculator* calculator;

};

#endif
