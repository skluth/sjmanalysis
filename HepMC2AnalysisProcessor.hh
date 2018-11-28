#ifndef HEPMC2ANALYSISPROCESSOR_HH
#define HEPMC2ANALYSISPROCESSOR_HH

#include "HepMC2Reader.hh"
#include <string>


class HepMC2AnalysisProcessor {

public:
  
  HepMC2AnalysisProcessor( const std::string & );

  void runAnalysis( Int_t maxevt=1000 );

private:

  HepMC2Reader reader;

};


#endif
