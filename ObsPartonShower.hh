#ifndef OBSPARTONSHOWER_HH
#define OBSPARTONSHOWER_HH

#include "Rtypes.h"
#include "ObsDifferential.hh"
#include "Analysis.hh"

#include <vector>
using std::vector;

class NtupleReader;
class TLorentzVector;
class DifferentialDataStructure;
class FilledObservable;

class ObsPartonShower : public ObsDifferential {

public:

  ObsPartonShower( const vector<Double_t>&, 
		   const vector<Double_t>&,
		   const vector<Double_t>&,
		   const vector<Analysis>&,
		   Double_t y34c=0.0045 );
  ~ObsPartonShower() {}
  virtual void addAnalyses( const vector<Analysis>& );
  virtual void fill( NtupleReader* ntr, const Analysis& variation );
  virtual vector<FilledObservable*> getFilledObservables() const;

private:

  vector<Double_t> getValues( NtupleReader* ntr, const string& reco );
  Double_t calcAngle( const TLorentzVector&, const TLorentzVector& );

  vector<Double_t> c202binedges;
  vector<Double_t> asbinedges;
  Double_t y34cut;
  map<string,DataStructure*> c202data;
  map<string,DataStructure*> asdata;

};

#endif
