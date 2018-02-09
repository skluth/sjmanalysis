#ifndef LEP2NTUPLEREADER_HH
#define LEP2NTUPLEREADER_HH

#include "NtupleReader.hh"

class LEP2NtupleReader : public NtupleReader {

  typedef std::pair<float,float> range;
  typedef std::map<std::string,range> rangemap;
  rangemap ecmsranges;

public:

  LEP2NtupleReader() {}
  LEP2NtupleReader( const char* filename, const char* ntid="h10", 
		    const bool lpr=true );
  virtual ~LEP2NtupleReader() {}

  virtual bool Preselection( const std::string& ecms );
  virtual bool Selection( const std::string& ecms );
  virtual const std::map<std::string,bool> getSelections( const std::string& ecms );

};

#endif
