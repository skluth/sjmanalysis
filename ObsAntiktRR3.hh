#ifndef OBSANTIKTRR3_HH
#define OBSANTIKTRR3_HH

#include "Observable.hh"
//#include "JetrateDataStructure.hh"
#include "Analysis.hh"

#include <vector>
#include <string>
//#include <map>
using std::vector;
using std::string;
//using std::map;

class NtupleReader;

class ObsAntiktRR3 : public Observable {

public:

  ObsAntiktRR3( const vector<Analysis>& variations, Double_t eminfrac=0.06 );
  ~ObsAntiktRR3() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation );

private:

  Double_t EminFraction;
  vector<Double_t> Rvalues;

};


#endif
