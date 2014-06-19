#ifndef OBSERVABLEFACTORY_HH
#define OBSERVABLEFACTORY_HH

#include "Rtypes.h"
#include "Analysis.hh"
#include <vector>
using std::vector;
#include <string>
using std::string;

class Observable;

class ObservableFactory {

public:

  ObservableFactory(); 
  ~ObservableFactory() {}
  
  vector<Observable*> createObservables( const vector<string>& obsnames,
					 const vector<Analysis>& analyses );

private:

  vector<Double_t> thrustbins;
  vector<Double_t> yNMbins;
  vector<Double_t> eminFraction;
  vector<Double_t> Rvalues;

};

#endif