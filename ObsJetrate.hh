#ifndef OBSJETRATE_HH
#define OBSJETRATE_HH

#include "Observable.hh"

#include "Rtypes.h"
#include "Analysis.hh"
#include <vector>
#include <string>
#include <map>

class NtupleReader;
class JetrateDataStructure;
class FilledObservable;
class JetrateCalculator;

class ObsJetrate : public Observable {

public:

  ObsJetrate( std::string name, 
	      const std::vector<Double_t> & pts,
	      const std::vector<Analysis> & variations,
	      const JetrateCalculator* calc,
	      const bool lprint=true );
  ObsJetrate();
  virtual ~ObsJetrate();
  virtual void fill( NtupleReader* ntr, const Analysis & variation );
  virtual std::vector<FilledObservable*> getFilledObservables() const;
  virtual void print() const;

private:

  virtual void addAnalysis( const Analysis & );
  void printDatastructures( const std::map<std::string,JetrateDataStructure*> & ) const;

  std::vector<Double_t> points;
  std::map<std::string,JetrateDataStructure*> jetrates2;
  std::map<std::string,JetrateDataStructure*> jetrates3;
  std::map<std::string,JetrateDataStructure*> jetrates4;
  std::map<std::string,JetrateDataStructure*> jetrates5;
  std::map<std::string,JetrateDataStructure*> jetrates6;
  const JetrateCalculator* calculator;

};

#endif
