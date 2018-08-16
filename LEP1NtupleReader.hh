#ifndef LEP1NTUPLEREADERLEP1_HH
#define LEP1NTUPLEREADERLEP1_HH

#include "LEPNtupleReader.hh"

class LEP1NtupleReader : public LEPNtupleReader {

public:

  LEP1NtupleReader() {}
  LEP1NtupleReader( const char* filename, const char* ntid="h10", 
		    const bool lpr=true );
  virtual ~LEP1NtupleReader() {}

  virtual bool Preselection( const std::string& ecms );
  virtual bool Selection( const std::string& ecms );
  virtual const std::map<std::string,bool> getSelections( const std::string& ecms );

};

#endif
