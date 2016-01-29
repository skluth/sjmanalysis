#ifndef OBSPARTONSHOWER_HH
#define OBSPARTONSHOWER_HH

#include "Rtypes.h"
#include "ObsDifferential.hh"
#include "Analysis.hh"

#include <vector>
using std::vector;
#include <iostream>

class NtupleReader;
class TLorentzVector;
class DifferentialDataStructure;
class FilledObservable;

class ObsPartonShower : public Observable {

public:

  ObsPartonShower( const vector<Double_t>&, 
		   const vector<Double_t>&,
		   const vector<Double_t>&,
		   const vector<Double_t>&,
		   const vector<Analysis>&,
		   Double_t y34c=0.0045, Double_t y34y23c=0.5 );
  virtual ~ObsPartonShower() {}
  virtual void fill( NtupleReader*, const Analysis& );
  virtual vector<FilledObservable*> getFilledObservables() const;

private:

  virtual void addAnalysis( const Analysis& );

  vector<Double_t> getValues( NtupleReader*, const string& );
  Double_t calcAngle( const TLorentzVector&, const TLorentzVector& );
  Double_t EnergyCorrelator( const vector<TLorentzVector>&, 
			     Int_t N=2, Double_t beta=0.2 );

  vector<Double_t> a14binedges;
  vector<Double_t> c202binedges;
  vector<Double_t> asbinedges;
  vector<Double_t> mrbinedges;
  Double_t y34cut;
  Double_t y34y23cut;
  map<string,DifferentialDataStructure*> a14data;
  map<string,DifferentialDataStructure*> c202data;
  map<string,DifferentialDataStructure*> asdata;
  map<string,DifferentialDataStructure*> mrdata;
  map<string,MatrixDataStructure*> a14matrices;
  map<string,MatrixDataStructure*> c202matrices;
  map<string,MatrixDataStructure*> asmatrices;
  map<string,MatrixDataStructure*> mrmatrices;

};

#endif
