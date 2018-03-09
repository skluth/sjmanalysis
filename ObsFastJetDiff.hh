#ifndef OBSFASTJETDIFF_HH
#define OBSFASTJETDIFF_HH

#include "Rtypes.h"
#include "ObsDifferential.hh"
#include "Analysis.hh"

#include <vector>
#include <string>
#include <map>

class NtupleReader;
class FilledObservable;
class DifferentialDataStructure;

class ObsFastJetDiff : public Observable {

public:

  ObsFastJetDiff( const std::string & name, const std::string & algo,
		  const std::vector<Double_t> & ynmbins, 
		  const std::vector<Analysis> & variations,
		  const bool lprint=true );
  ~ObsFastJetDiff();
  virtual void fill( NtupleReader* ntr, const Analysis & variation );
  virtual std::vector<FilledObservable*> getFilledObservables() const;

private:

  virtual void addAnalysis( const Analysis & );

  std::string Algorithm;
  std::vector<Double_t> binedges;
  std::map<std::string,DifferentialDataStructure*> ymerge23;
  std::map<std::string,DifferentialDataStructure*> ymerge34;
  std::map<std::string,DifferentialDataStructure*> ymerge45;
  std::map<std::string,DifferentialDataStructure*> ymerge56;

};

#endif
