#ifndef OBSANTIKTEMINR3_HH
#define OBSANTIKTEMINR3_HH

#include "Observable.hh"
#include "Analysis.hh"

#include <vector>
#include <string>
using std::vector;
using std::string;

class NtupleReader;

class ObsAntiktEminR3 : public Observable {

public:

  ObsAntiktEminR3( const vector<Analysis>& variations, Double_t R=0.7 );
  ~ObsAntiktEminR3() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation );

private:

  Double_t Remin;
  vector<Double_t> eminFraction;

};


#endif
