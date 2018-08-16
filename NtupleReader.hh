#ifndef NTUPLEREADER_HH
#define NTUPLEREADER_HH

#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TString.h"
#include <vector>
#include <string>
#include <map>

class NtupleReader {

public:

  NtupleReader() {}
  virtual ~NtupleReader() {}

  virtual Int_t GetNumberEntries() = 0;
  virtual bool GetEvent( Int_t ievnt ) = 0;
 
  virtual const std::vector<TLorentzVector> GetLorentzVectors( const std::string & opt ) = 0;
  
  virtual Double_t Evis( const std::vector<TLorentzVector>& v ) const = 0;
  
  virtual const std::map<std::string,bool> getSelections( const std::string& ) = 0;

  virtual bool MCNonRad() = 0;
  virtual bool isMC() = 0;

  virtual Float_t abscostt() = 0;
  // virtual Double_t getYmerge( const TString& algorithm, const TString& reco, Int_t njet ) = 0;
  // virtual Double_t getThrust( const TString& reco ) = 0;
  
};

#endif
