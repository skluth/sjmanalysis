#ifndef OBSEEC_HH
#define OBSEEC_HH

#include "Rtypes.h"
#include "Observable.hh"
#include "Analysis.hh"

#include <iostream>
#include <vector>
#include <string>
#include <map>

class NtupleReader;
class DifferentialDataStructure;
class TLorentzVector;

class ObsEEC : public Observable {

public:

  ObsEEC( const std::string & name,
	  const std::vector<Double_t> & bins,
	  const std::vector<Analysis> & variations,
	  const bool scOpt=true,
	  const bool lprint=true );
  virtual ~ObsEEC();
  virtual std::vector<FilledObservable*> getFilledObservables() const;
  virtual void fill( NtupleReader* ntr, const Analysis & variation );

protected:

  virtual void calcWeight( const TLorentzVector&, const TLorentzVector&,
                           Double_t& , Double_t& );

  bool selfCorrelation;

private:

  virtual void addAnalysis( const Analysis & );

  std::vector<Double_t> binedges;
  std::map<std::string,DifferentialDataStructure*> data;
  std::map<std::string,DifferentialDataStructure*> datab;
  std::map<std::string,DifferentialDataStructure*> dataudsc;

  Int_t nevents;
  Int_t neventsb;
  Int_t neventsudsc;
  
};

#endif
