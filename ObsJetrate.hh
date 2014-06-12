#ifndef OBSJETRATE_HH
#define OBSJETRATE_HH

#include "Observable.hh"
#include "Rtypes.h"
#include "Analysis.hh"

#include <vector>
#include <string>
#include <map>
using std::vector;
using std::string;
using std::map;

class NtupleReader;
class DataStructure;
//class JetrateDataStructure;

class ObsJetrate : public Observable {

public:

  ObsJetrate( string name );
  ObsJetrate() {}
  ~ObsJetrate() {}
  void addAnalyses( const vector<Analysis>& variations, const vector<Double_t>& points );
  virtual void fill( NtupleReader* ntr, const Analysis& variation ) = 0;
  void finalise();
  void print();
  // virtual map<string,DataStructure*> getData();
  //virtual DataStructure getDataStructure( const Analysis& );
  //virtual void setDataStructure( DataStructure*, const Analysis& );

protected:

  //  map<string,DataStructure*> datastructures;

};


#endif
