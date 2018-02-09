#ifndef OBSERVABLEFACTORY_HH
#define OBSERVABLEFACTORY_HH

#include "Rtypes.h"
#include "Analysis.hh"
#include "SjmConfigParser.hh"
#include <vector>
using std::vector;
#include <string>
using std::string;

class Observable;

class ObservableFactory {

  SjmConfigParser sjmConfigs;

  bool nameIs( const string&, const string& );

public:

  ObservableFactory( const SjmConfigParser& ); 
  ~ObservableFactory() {}
  
  vector<Observable*> createObservables( const vector<string>& obsnames,
					 const vector<Analysis>& analyses );

};

#endif
