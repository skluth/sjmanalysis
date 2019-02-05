#ifndef HEPMCROOTREADER_HH
#define HEPMCROOTREADER_HH

#include "NtupleReader.hh"

#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TString.h"
#include <vector>
#include <string>
#include <map>
#include <utility>

class TFile;
class TTree;
namespace HepMC {
  class GenEventData;
  class GenRunInfoData;
  class GenParticleData;
  class GenVertexData;
}

class HepMCRootReader : public NtupleReader {

public:

  HepMCRootReader( const std::string & filename,
		   const std::string & treename="h10",
		   const std::string & branchname="h10" );
  virtual ~HepMCRootReader();

  virtual Int_t GetNumberEntries();
  virtual bool GetEvent( Int_t ievnt );
  virtual bool GetNextEvent( Int_t maxevt=0 );
 
  virtual const std::vector<TLorentzVector> GetLorentzVectors( const std::string & opt );
    
  virtual const std::map<std::string,bool> getSelections( const std::string & );

  virtual bool MCNonRad();
  virtual bool isMC();
  virtual void printParticlesVertices();

  typedef HepMC::GenParticleData Particle;
  typedef HepMC::GenVertexData Vertex;
  typedef std::vector<const Particle*> ParticleVector;
  typedef std::map<const Vertex*,ParticleVector> VertexParticlesMap;
  typedef std::map<const Particle*,const Vertex*> ParticleVertexMap;
  
private:

  void getHadron();
  void getParton();
  void getIsr();

  Double_t event_scale();
  ParticleVector getBeamParticles();
  void findISRphotons();
  
  // Variables for reading
  Int_t nevents;
  TFile* hepmcrootfile;
  TTree* hepmctree;
  HepMC::GenEventData* eventData;
  HepMC::GenRunInfoData* runInfoData;
  VertexParticlesMap incomingParticles;
  VertexParticlesMap outgoingParticles;
  ParticleVertexMap beginVertexMap;
  ParticleVertexMap endVertexMap;
  
  // Variables for returning objects:
  std::vector<TLorentzVector> vtlv;
  std::map<std::string,std::vector<TLorentzVector>> vtlvCache;
  std::map<std::string,Bool_t> cacheIsValid;
  
  ParticleVector ISRphotons;

};


//Double_t Evis( const std::vector<TLorentzVector> & v );


#endif
