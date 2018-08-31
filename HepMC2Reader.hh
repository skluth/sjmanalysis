#ifndef HEPMC2READER_HH
#define HEPMC2READER_HH

#include "NtupleReader.hh"
#include "TLorentzVector.h"
#include "TString.h"
#include <vector>
#include <string>
#include <map>
#include <iostream>

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"


class HepMC2Reader : public NtupleReader {

public:

  HepMC2Reader( const std::string& filename );
  virtual ~HepMC2Reader() {}

  virtual Int_t GetNumberEntries() { return 0; }
  virtual bool GetEvent( Int_t ievnt=0 );
 
  virtual const std::vector<TLorentzVector> GetLorentzVectors( const std::string & opt );
  
  virtual const std::map<std::string,bool> getSelections( const std::string& ) {
    return std::map<std::string,bool>();
  }

  virtual bool MCNonRad() { return true; }
  virtual bool isMC() { return true; }

private:

  std::ifstream input;
  HepMC::GenEvent event;
  
};

#endif
