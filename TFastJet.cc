
#include "TFastJet.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/SISConePlugin.hh"
#include "fastjet/PxConePlugin.hh"
#include "fastjet/SISConeSphericalPlugin.hh"
#include "fastjet/JadePlugin.hh"

#include "TParticle.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <map>
using std::map;
#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;
#include <stdexcept>

using std::vector;

// Declare Recombiner class for E0 scheme e+e- algorithms:
class EEE0Recombiner: public fastjet::JetDefinition::Recombiner {
public:
  EEE0Recombiner() {}
  ~EEE0Recombiner() {}
  virtual string description() const;
  virtual void recombine( const fastjet::PseudoJet& pa,
			  const fastjet::PseudoJet& pb, 
			  fastjet::PseudoJet& pab ) const;
};

// Wrapper/adapter to fastjet classes:
TFastJet::TFastJet( const vector<TLorentzVector> & vtl,
		    const char* jetalg,
		    const double R,
		    const double Emin ) : clusseq(0) {

  // Map of supported algorithms to fastjet enum:
  map<string,fastjet::JetAlgorithm> jamap {
    { "kt", fastjet::kt_algorithm },
    { "cambridge", fastjet::cambridge_algorithm },
    { "antikt", fastjet::antikt_algorithm },
    { "genkt", fastjet::genkt_algorithm },
    { "siscone", fastjet::plugin_algorithm },
    { "eesiscone", fastjet::plugin_algorithm },
    { "eekt", fastjet::ee_kt_algorithm },
    { "jade", fastjet::plugin_algorithm },
    { "eeantikt", fastjet::ee_genkt_algorithm },
    { "pxcone", fastjet::plugin_algorithm }
  };
  string jetalgString( jetalg );
  if( jamap.find( jetalgString ) == jamap.end() ) {
    string txt= "TFastJet::TFastJet: wrong algorithm: " + jetalgString;
    throw std::runtime_error( txt );
  }

  // Setup fastjet components, there is no mem leak for Recombiner and Plugin
  // objects since ownership is passed to fastjet::JetDefinition object:
  fastjet::JetAlgorithm ja= jamap[jetalgString];
  fastjet::JetDefinition jetdef;
  fastjet::JetDefinition::Recombiner* recombiner= 0;
  fastjet::JetDefinition::Plugin* plugin= 0;
  if( ja == fastjet::plugin_algorithm ) {
    if( jetalgString == "siscone" ) {
      plugin= new fastjet::SISConePlugin( R, 0.75 );
    }
    else if( jetalgString == "eesiscone" ) {
      plugin= new fastjet::SISConeSphericalPlugin( R, 0.75 );
    }
    else if( jetalgString == "jade" ) {
      plugin= new fastjet::JadePlugin();
      recombiner= new EEE0Recombiner();
    }
    else if( jetalgString == "pxcone" ) {
      plugin= new fastjet::PxConePlugin( R, Emin, 0.75, true, 1 );
    }
    else {
      string txt= "TFastJet::TFastJet: wrong plugin: " + jetalgString;
      throw std::runtime_error( txt );
    }
    jetdef= fastjet::JetDefinition( plugin );
    jetdef.delete_plugin_when_unused();
    if( recombiner ) {
      jetdef.set_recombiner( recombiner );
      jetdef.delete_recombiner_when_unused();
    }
  }
  else if( ja == jamap["eekt"] ) {
    jetdef= fastjet::JetDefinition( ja );
  }
  else if( ja == jamap["eeantikt"] ) {
    jetdef= fastjet::JetDefinition( ja, R, -1.0 );
  }
  else {
    jetdef= fastjet::JetDefinition( ja, R );
  }

  // Run jet algorithm:
  clusseq= new fastjet::ClusterSequence( vtl, jetdef );

}

TFastJet::~TFastJet() {
  if( clusseq ) delete clusseq;
}

// All jets with p_t > p_t,min:
const vector<TLorentzVector> TFastJet::inclusive_jets( const double ptmin ) {
  vector<fastjet::PseudoJet> incljets= clusseq->inclusive_jets( ptmin );
  vector<fastjet::PseudoJet> pjets= sorted_by_pt( incljets );
  return copyPseudoJetsToLorentzVectors( pjets );
}

// All jets with E > Emin for e+e-:
const vector<TLorentzVector> TFastJet::inclusive_eejets( const double Emin ) {
  vector<fastjet::PseudoJet> incljets= clusseq->inclusive_jets( 0.0 );
  vector<fastjet::PseudoJet> pjets= sorted_by_E( incljets );
  return copyPseudoJetsToLorentzVectors( pjets, Emin );
}

// Fixed number of jets (e+e-):
const vector<TLorentzVector> TFastJet::exclusive_eejets( const int njets ) {
  vector<fastjet::PseudoJet> excljets= clusseq->exclusive_jets( njets );
  vector<fastjet::PseudoJet> pjets= sorted_by_E( excljets );
  return copyPseudoJetsToLorentzVectors( pjets );
}

// Same with index maps of jets to input 4-vectors:
const vector<TLorentzVector>
TFastJet::inclusive_eejets( vector< vector<int> > & vindx,
			    const double Emin ) {
  vector<fastjet::PseudoJet> incljets= clusseq->inclusive_jets( 0.0 );
  vector<fastjet::PseudoJet> pjets= sorted_by_E( incljets );
  vindx= constituents( pjets );
  return copyPseudoJetsToLorentzVectors( pjets, Emin );
}
const vector<TLorentzVector>
TFastJet::inclusive_jets( vector< vector<int> > & vindx,
			  const double ptmin ) {
  vector<fastjet::PseudoJet> incljets= clusseq->inclusive_jets( ptmin );
  vector<fastjet::PseudoJet> pjets= sorted_by_pt( incljets );
  vindx= constituents( pjets );
  return copyPseudoJetsToLorentzVectors( pjets );
}
const vector<TLorentzVector>
TFastJet::exclusive_eejets( vector< vector<int> > & vindx,
			    const int njets ) {
  vector<fastjet::PseudoJet> excljets= clusseq->exclusive_jets( njets );
  vector<fastjet::PseudoJet> pjets= sorted_by_E( excljets );
  vindx= constituents( pjets );
  return copyPseudoJetsToLorentzVectors( pjets );
}

// Internal helper:
const vector<TLorentzVector>
TFastJet::copyPseudoJetsToLorentzVectors( const vector<fastjet::PseudoJet> & pjets,
					  const double Emin ) {
  vector<TLorentzVector> jetstlv;
  jetstlv.reserve( pjets.size() );
  for( const fastjet::PseudoJet & pj : pjets ) {
    if( pj.E() > Emin ) {
      TLorentzVector tlv( pj.px(), pj.py(), pj.pz(), pj.E() );
      jetstlv.push_back( tlv );
    }
    else {
      break;
    }
  }
  return jetstlv;
}

// Jet constituents index map into array of input pseudojets:
vector< vector<int> > TFastJet::constituents( const vector<fastjet::PseudoJet> & pjets ) {
  vector<int> p2j= clusseq->particle_jet_indices( pjets );
  vector< vector<int> > constituentsMaps;
  for( size_t ijet= 0; ijet < pjets.size(); ijet++ ) {
    vector<int> constituentsMap;
    for( size_t i= 0; i < p2j.size(); i++ ) {
      if( p2j[i] == int(ijet) ) {
	constituentsMap.push_back( i );
      }
    }
    constituentsMaps.push_back( constituentsMap );
  }
  return constituentsMaps;
}

// Ymerge for transition n+1 -> n jets (e+e-):
double TFastJet::ymerge( int njets ) {
  return clusseq->exclusive_ymerge_max( njets );
}

// Number of jets when all ymerge > ycut (e+e-):
int TFastJet::njets( double ycut ) {
  return clusseq->n_exclusive_jets_ycut( ycut );
}

// Visible energy a la fastjet:
double TFastJet::Evis() {
  return clusseq->Q();
}

// class EEE0Recombiner implementation to study remobination effects with Jade:
string EEE0Recombiner::description() const {
  return "E0 scheme for EE";
}

void EEE0Recombiner::recombine( const fastjet::PseudoJet & pa,
				const fastjet::PseudoJet & pb,
				fastjet::PseudoJet & pab ) const {
  pab.reset( pa.px() + pb.px(), pa.py() + pb.py(),
	     pa.pz() + pb.pz(), pa.E() + pb.E() );
  // double rescale= pab.E()/TMath::Sqrt( pab.perp2() + pab.pz()*pab.pz() );
  // pab.reset_momentum( rescale*pab.px(), rescale*pab.py(), rescale*pab.pz(), pab.E() );
  return;
}

