#ifndef OBSGROOMED_HH
#define OBSGROOMED_HH

#include "Rtypes.h"
#include "ObsDifferential.hh"
#include "Analysis.hh"

#include <vector>
#include <map>
#include <functional>

class string;
class NtupleReader;
class TLorentzVector;
class DifferentialDataStructure;
class FilledObservable;

class ObsGroomed : public Observable {

public:

  ObsGroomed( const std::vector<Double_t> &,
	      const std::vector<Double_t> &,
	      const std::vector<Analysis> &,
	      bool lprint=true );
  virtual ~ObsGroomed() {}
  virtual void fill( NtupleReader*, const Analysis & );
  virtual std::vector<FilledObservable*> getFilledObservables() const;

private:
  
  virtual void addAnalysis( const Analysis & );
  const std::map<std::string,Double_t> getValues( NtupleReader*, const std::string & ) const;
  typedef std::function< void ( const std::string &, const std::string & ) > LoopFunc;
  void loop( const LoopFunc & ) const;
  
  std::map<std::string,std::vector<Double_t>> binedges;
  std::vector<Double_t> betaValues;
  std::vector<Double_t> zcutValues;

  typedef std::map<std::string,DifferentialDataStructure*> DdsMap;
  typedef std::map<std::string,MatrixDataStructure*> MdsMap;
  std::map<std::string,DdsMap> grData;
  std::map<std::string,MdsMap> grMatrices;

};

#endif
