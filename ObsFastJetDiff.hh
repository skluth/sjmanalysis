#ifndef OBSFASTJETDIFF_HH
#define OBSFASTJETDIFF_HH

#include "Rtypes.h"
#include "ObsDifferential.hh"
#include "Analysis.hh"

#include <vector>
using std::string;
#include <string>
using std::vector;
#include <map>
using std::map;

class NtupleReader;
class FilledObservable;
class DataStructure;

class ObsFastJetDiff : public Observable {

public:

  ObsFastJetDiff( const string& name, const string& algo,
		  const vector<Double_t>& ynmbins, 
		  const vector<Analysis>& variations,
		  const bool lprint=true );
  ~ObsFastJetDiff();
  virtual void fill( NtupleReader* ntr, const Analysis& variation );
  virtual vector<FilledObservable*> getFilledObservables() const;
  virtual bool containsAnalysis( const Analysis& );  
  virtual void addAnalyses( const vector<Analysis>& );

private:

  string Algorithm;
  vector<Double_t> binedges;
  map<string,DifferentialDataStructure*> ymerge23;
  map<string,DifferentialDataStructure*> ymerge34;
  map<string,DifferentialDataStructure*> ymerge45;
  map<string,DifferentialDataStructure*> ymerge56;

};

#endif
