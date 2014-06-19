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

class ObsFastJetDiff : public ObsDifferential {

public:

  ObsFastJetDiff( const string& name, const string& algo,
		  const vector<Double_t>& bins, 
		  const vector<Analysis>& variations );
  ~ObsFastJetDiff();
  virtual void addAnalyses( const vector<Analysis>& variations );
  virtual void fill( NtupleReader* ntr, const Analysis& variation );
  virtual vector<FilledObservable*> getFilledObservables() const;

private:
  
  string Algorithm;
  map<string,DataStructure*> ymerge23;
  map<string,DataStructure*> ymerge34;
  map<string,DataStructure*> ymerge45;
  map<string,DataStructure*> ymerge56;

};

#endif
