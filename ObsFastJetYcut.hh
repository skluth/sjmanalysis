#ifndef OBSFASTJETYCUT_HH
#define OBSFASTJETYCUT_HH

#include "ObsJetrate.hh"
#include "Analysis.hh"

#include <vector>
using std::vector;
#include <string>
using std::string;

class NtupleReader;

class ObsFastJetYcut : public ObsJetrate {

public:

  ObsFastJetYcut( const string& name, const string& algo,
		  const vector<Double_t>& ycuts, 
		  const vector<Analysis>& variations );
  ~ObsFastJetYcut() {}
  virtual void fill( NtupleReader* ntr, const Analysis& variation );
  
private:
  
  string Algorithm;
  vector<Double_t> YcutPoints;

};

#endif
