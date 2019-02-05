#ifndef HEPMC2READER_HH
#define HEPMC2READER_HH

#include "NtupleReader.hh"
#include "TLorentzVector.h"
#include <vector>
#include <string>
#include <map>
#include <iostream>

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/IO_GenEvent.h"


class HepMC2Reader : public NtupleReader {

public:

  HepMC2Reader( const std::string& filename );
  virtual ~HepMC2Reader();

  virtual Int_t GetNumberEntries() { return 0; }
  virtual bool GetEvent( Int_t ievent=0 );
  virtual bool GetNextEvent( Int_t maxevt=0 );
 
  virtual const std::vector<TLorentzVector> GetLorentzVectors( const std::string & opt );
  
  virtual const std::map<std::string,bool> getSelections( const std::string & ) {
    return std::map<std::string,bool>();
  }

  virtual bool MCNonRad();
  virtual bool isMC() { return true; }

  virtual void printParticlesVertices();

  
private:
  
  void findISRphotons();
  void getHadron();
  void getParton();
  void getIsr();

  std::ifstream input;
  Int_t nevent;
  HepMC::GenEvent event;
  std::vector<const HepMC::GenParticle*> ISRphotons;

  std::vector<TLorentzVector> vtlv;
  std::map<std::string,std::vector<TLorentzVector>> vtlvCache;
  std::map<std::string,Bool_t> cacheIsValid;

};

#endif
