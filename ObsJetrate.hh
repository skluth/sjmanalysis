#ifndef OBSJETRATE_HH
#define OBSJETRATE_HH

#include "Observable.hh"

#include "Rtypes.h"
#include "Analysis.hh"
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <map>
using std::map;

class NtupleReader;
class JetrateDataStructure;
class FilledObservable;
class JetrateCalculator;

class ObsJetrate : public Observable {

public:

  ObsJetrate( string, 
	      const vector<Double_t>&,
	      const vector<Analysis>& variations,
	      const JetrateCalculator* calc,
	      const bool lprint=true );
  ObsJetrate();
  virtual ~ObsJetrate();
  virtual void fill( NtupleReader* ntr, const Analysis& variation );
  virtual vector<FilledObservable*> getFilledObservables() const;
  virtual void print() const;

private:

  virtual void addAnalysis( const Analysis& );
  void printDatastructures( const map<string,JetrateDataStructure*>& ) const;

  vector<Double_t> points;
  map<string,JetrateDataStructure*> jetrates2;
  map<string,JetrateDataStructure*> jetrates3;
  map<string,JetrateDataStructure*> jetrates4;
  map<string,JetrateDataStructure*> jetrates5;
  map<string,JetrateDataStructure*> jetrates6;
  const JetrateCalculator* calculator;

};

#endif
