#ifndef OBSERVABLEFACTORY_HH
#define OBSERVABLEFACTORY_HH

#include "Rtypes.h"
#include "Analysis.hh"
#include "SjmConfigParser.hh"
#include <vector>
#include <string>

class Observable;

class ObservableFactory {

  SjmConfigParser sjmConfigs;

  bool nameIs( const std::string &, const std::string & );

public:

  ObservableFactory( const SjmConfigParser & ); 
  ~ObservableFactory() {}
  
  std::vector<Observable*>
  createObservables( const std::vector<std::string> & obsnames,
		     const std::vector<Analysis> & analyses );

};

#endif
