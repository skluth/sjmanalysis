
#include "TFastJet.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/SISConePlugin.hh"

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
using std::logic_error;

// TFastJet::TFastJet( const vector<TParticle>& vtp ) : clusseq(0),
// 						     plugin(0),
// 						     pjets(0) {
//   pjets= new vector<fastjet::PseudoJet>();
//   vector<fastjet::PseudoJet> particles;
//   for( UInt_t i= 0; i < vtp.size(); i++ ) {
//     const TParticle& tp= vtp[i];
//     fastjet::PseudoJet pj= fastjet::PseudoJet( tp.Px(), tp.Py(), 
// 					       tp.Pz(), tp.Energy() );
//     pj.set_user_index( i );
//     particles.push_back( pj );
//   }
//   double R= 0.4;
//   fastjet::JetDefinition jetdef( fastjet::antikt_algorithm, R );
//   clusseq= new fastjet::ClusterSequence( particles, jetdef );  
//   return;
// }

TFastJet::TFastJet( const vector<TLorentzVector>& vtl, 
		    const char* jetalg, 
		    const double& R, 
		    const vector<int>* vindx ) : clusseq(0),
						 plugin(0),
						 pjets(0) {
  // Hide map from rootcint
  static map<string,fastjet::JetAlgorithm> jamap;
  if( jamap.size() == 0 ) {
    jamap["kt"]= fastjet::kt_algorithm;
    jamap["cambridge"]= fastjet::cambridge_algorithm;
    jamap["antikt"]= fastjet::antikt_algorithm;
    jamap["genkt"]= fastjet::genkt_algorithm;
    jamap["siscone"]= fastjet::plugin_algorithm;
    jamap["eesiscone"]= fastjet::plugin_algorithm;
    jamap["eekt"]= fastjet::ee_kt_algorithm;
    jamap["jade"]= fastjet::plugin_algorithm;
    jamap["eeantikt"]= fastjet::ee_genkt_algorithm;
  }
  // Check input:
  string jetalgString( jetalg );
  if( jamap.find( jetalgString ) == jamap.end() ) {
    string txt= "TFastJet::TFastJet: wrong algorithm: " + jetalgString;
    throw logic_error( txt );
  }
  // Setup fastjet:
  fastjet::JetAlgorithm ja= jamap[jetalgString];
  fastjet::JetDefinition jetdef;
  fastjet::JetDefinition::Recombiner* recombiner= 0;
  //  if( ja == fastjet::plugin_algorithm ) {
  if( ja == 99 ) {
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
    else {
      string txt= "TFastJet::TFastJet: wrong plugin: " + jetalgString;
      throw logic_error( txt );
    }
    jetdef= fastjet::JetDefinition( plugin );
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
  // Prepare input 4-vectors:
  vector<fastjet::PseudoJet> particles;
  for( UInt_t i= 0; i < vtl.size(); i++ ) {
    fastjet::PseudoJet pj( vtl[i] );
    if( vindx == 0 ) pj.set_user_index( i );
    else pj.set_user_index( vindx->at(i) );
    particles.push_back( pj );
  }
  // Run jet algorithm:
  clusseq= new fastjet::ClusterSequence( particles, jetdef );
  // For receiving output from fastjet:
  pjets= new vector<fastjet::PseudoJet>();
}

TFastJet::~TFastJet() {
  if( plugin ) delete plugin;
  if( clusseq ) delete clusseq;
  if( pjets ) delete pjets;
}

// All jets with p_t > p_t,min (or E < Emin for e+e-):
const vector<TLorentzVector>& TFastJet::inclusive_jets( const double& ptmin ) {
  vector<fastjet::PseudoJet> incljets= clusseq->inclusive_jets( ptmin );
  *pjets= sorted_by_pt( incljets );
  return copyPseudoJetsToLorentzVectors();
}

// Fixed number of jets (e+e-):
const vector<TLorentzVector>& TFastJet::exclusive_jets( const int njets ) {
  vector<fastjet::PseudoJet> excljets= clusseq->exclusive_jets( njets );
  *pjets= sorted_by_E( excljets );
  return copyPseudoJetsToLorentzVectors();
}

// Internal helper:
const vector<TLorentzVector>& TFastJet::copyPseudoJetsToLorentzVectors() {
  static vector<TLorentzVector> jetstlv;
  jetstlv.clear();
  jetstlv.reserve( pjets->size() );
  for( UInt_t i= 0; i < pjets->size(); i++ ) {
    fastjet::PseudoJet pj= (*pjets)[i];
    TLorentzVector tlv( pj.px(), pj.py(), pj.pz(), pj.E() );
    jetstlv.push_back( tlv );
  }
  return jetstlv;
}

// Jet constituents index map into array of input 4-vectors:
vector< vector<int> >& TFastJet::constituents() {
  vector< vector<int> >* cnstmap= new vector< vector<int> >();
  fastjet::PseudoJet pj;
  for( UInt_t i= 0; i < pjets->size(); i++ ) {
    pj= (*pjets)[i];
    vector<int> vindx;
    vector<fastjet::PseudoJet> cnst= clusseq->constituents( pj );
    for( UInt_t j= 0; j < cnst.size(); j++ ) {
      vindx.push_back( cnst[j].user_index() );
    }
    cnstmap->push_back( vindx );
  }
  return *cnstmap;
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

// class EEE0Recombiner to study remobination effects with Jade:

string EEE0Recombiner::description() const {
  return "E0 scheme for EE"; 
}

void EEE0Recombiner::recombine( const fastjet::PseudoJet& pa,
				const fastjet::PseudoJet& pb, 
				fastjet::PseudoJet& pab ) const {
  pab.reset( pa.px() + pb.px(), pa.py() + pb.py(),
	     pa.pz() + pb.pz(), pa.E() + pb.E() );
  // double rescale= pab.E()/TMath::Sqrt( pab.perp2() + pab.pz()*pab.pz() );
  // pab.reset_momentum( rescale*pab.px(), rescale*pab.py(), rescale*pab.pz(), pab.E() );
  return;
}

