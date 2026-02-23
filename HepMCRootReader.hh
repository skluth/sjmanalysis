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
namespace HepMC3 {
  class GenEventData;
  class GenRunInfoData;
  class GenParticleData;
  class GenVertexData;
}

class HepMCRootReader : public NtupleReader {

public:

  HepMCRootReader( const std::string & filename,
		   const std::string & treename,
		   const std::string & branchname,
		   const Double_t mcnonradcutin=1.0 );
  virtual ~HepMCRootReader();

  virtual Int_t GetNumberEntries();
  virtual bool GetEvent( Int_t ievnt );
  virtual bool GetNextEvent( Int_t maxevt=0 );
 
  virtual const std::vector<TLorentzVector> GetLorentzVectors( const std::string & opt );
    
  virtual const std::map<std::string,bool> getSelections( const std::string & );

  virtual bool MCNonRad();
  virtual bool isMC();
  virtual Int_t getPrimaryFlavour();
  virtual void printParticlesVertices();

  typedef HepMC3::GenParticleData Particle;
  typedef HepMC3::GenVertexData Vertex;
  typedef std::vector<const Particle*> ParticleVector;
  typedef std::map<const Vertex*,ParticleVector> VertexParticlesMap;
  typedef std::map<const Particle*,const Vertex*> ParticleVertexMap;
  
private:

  void getHadron();
  void getParton();
  void getIsr();
  
  std::vector<size_t> getIndices( const ParticleVector& ) const;
  ParticleVector getParents( const Particle* particle ) const;
  ParticleVector getDaughters( const Particle* particle ) const;

  Double_t event_scale();
  ParticleVector getBeamParticles();
  void findISRphotons();
  
  // Variables for reading
  Int_t nevents;
  TFile* hepmcrootfile;
  TTree* hepmctree;
  Double_t mcnonradcut;
  HepMC3::GenEventData* eventData;
  HepMC3::GenRunInfoData* runInfoData;
  VertexParticlesMap incomingParticlesMap;
  VertexParticlesMap outgoingParticlesMap;
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
