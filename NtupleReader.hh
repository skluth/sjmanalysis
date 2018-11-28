#ifndef NTUPLEREADER_HH
#define NTUPLEREADER_HH

#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TString.h"
#include <vector>
#include <string>
#include <map>
#include <algorithm>

class NtupleReader {

public:

  NtupleReader() {}
  virtual ~NtupleReader() {}

  virtual Int_t GetNumberEntries() = 0;
  virtual bool GetEvent( Int_t ievnt ) = 0;
  virtual bool GetNextEvent( Int_t maxevt=0 ) = 0;
 
  virtual const std::vector<TLorentzVector> GetLorentzVectors( const std::string & opt ) = 0;
    
  virtual const std::map<std::string,bool> getSelections( const std::string& ) = 0;

  virtual bool MCNonRad() = 0;
  virtual bool isMC() = 0;
  
};


Double_t Evis( const std::vector<TLorentzVector>& v );


#endif
