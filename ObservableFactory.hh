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

  Observable* createThrust( const vector<Analysis>& analyses );
  Observable* createFastJetYcut( const string& obsname, 
				 const string& algo,
				 const vector<Analysis>& analyses );
  Observable* createFastJetEmin( const string& obsname,
				 const string& algo,
				 const vector<Analysis>& analyses,
				 Double_t rvalue=0.7 );
  Observable* createFastJetR( const string& obsname,
			      const string& algo,
			      const vector<Analysis>& analyses,
			      Double_t eminfrac=0.06 );
  Int_t getNjetFromName( const string& name );

  vector<Double_t> yNMbins;

};

#endif
