#ifndef OBSPARTONSHOWER_HH
#define OBSPARTONSHOWER_HH

#include "Rtypes.h"
#include "ObsDifferential.hh"
#include "Analysis.hh"

#include <vector>
#include <map>
#include <iostream>

class NtupleReader;
class TLorentzVector;
class DifferentialDataStructure;
class FilledObservable;

class ObsPartonShower : public Observable {

public:

  ObsPartonShower( const std::vector<Double_t> &,
		   const std::vector<Double_t> &,
		   const std::vector<Double_t> &,
		   const std::vector<Double_t> &,
		   const std::vector<Analysis> &,
		   Double_t y34c=0.0045, Double_t y34y23c=0.5,
		   bool lprint=true );
  virtual ~ObsPartonShower() {}
  virtual void fill( NtupleReader*, const Analysis & );
  virtual std::vector<FilledObservable*> getFilledObservables() const;

private:

  virtual void addAnalysis( const Analysis & );

  std::vector<Double_t> getValues( NtupleReader*, const std::string & );
  Double_t calcAngle( const TLorentzVector &, const TLorentzVector & );
  Double_t EnergyCorrelator( const std::vector<TLorentzVector> &,
			     Int_t N=2, Double_t beta=0.2 );

  std::vector<Double_t> a14binedges;
  std::vector<Double_t> c202binedges;
  std::vector<Double_t> asbinedges;
  std::vector<Double_t> mrbinedges;
  Double_t y34cut;
  Double_t y34y23cut;
  std::map<std::string,DifferentialDataStructure*> a14data;
  std::map<std::string,DifferentialDataStructure*> c202data;
  std::map<std::string,DifferentialDataStructure*> asdata;
  std::map<std::string,DifferentialDataStructure*> mrdata;
  std::map<std::string,MatrixDataStructure*> a14matrices;
  std::map<std::string,MatrixDataStructure*> c202matrices;
  std::map<std::string,MatrixDataStructure*> asmatrices;
  std::map<std::string,MatrixDataStructure*> mrmatrices;

};

#endif
