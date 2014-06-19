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
class DataStructure;
class JetrateDataStructure;
class FilledObservable;

class ObsJetrate : public Observable {

public:

  ObsJetrate( string, const vector<Double_t>& );
  ObsJetrate() {}
  virtual ~ObsJetrate();
  virtual void addAnalyses( const vector<Analysis>& variations );
  virtual void fill( NtupleReader* ntr, const Analysis& variation ) = 0;
  virtual vector<FilledObservable*> getFilledObservables() const;
  virtual void print() const;
  virtual bool containsAnalysis( const Analysis& );

protected:

  void getAndFillJetrateDataStructures( const vector<Double_t>& NJets, 
   					const string& tag );

  vector<Double_t> points;

private:

  JetrateDataStructure* getJetrateDataStructure( DataStructure* ds ) const;
  void getAndFillJetrateDataStructure( const vector<Double_t>&, 
				       const string&,
				       const map<string,DataStructure*>& );
  void printDatastructures( const map<string,DataStructure*>& ) const;

  map<string,DataStructure*> jetrates2;
  map<string,DataStructure*> jetrates3;
  map<string,DataStructure*> jetrates4;
  map<string,DataStructure*> jetrates5;
  map<string,DataStructure*> jetrates6;

};

#endif
