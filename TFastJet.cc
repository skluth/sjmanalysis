
#include "TFastJet.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/SISConePlugin.hh"
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

TFastJet::TFastJet( const vector<TParticle>& vtp ) : clusseq(0),
						     plugin(0),
						     pjets(0) {
  pjets= new vector<fastjet::PseudoJet>();
  vector<fastjet::PseudoJet> particles;
  for( UInt_t i= 0; i < vtp.size(); i++ ) {
    const TParticle& tp= vtp[i];
    fastjet::PseudoJet pj= fastjet::PseudoJet( tp.Px(), tp.Py(), 
					       tp.Pz(), tp.Energy() );
    pj.set_user_index( i );
    particles.push_back( pj );
  }
  double R= 0.4;
  fastjet::JetDefinition jetdef( fastjet::antikt_algorithm, R );
  clusseq= new fastjet::ClusterSequence( particles, jetdef );  
  return;
}

TFastJet::TFastJet( const vector<TLorentzVector>& vtl, 
		    const char* jetalg, 
		    const double& R, 
		    const vector<int>* vindx ) : clusseq(0),
						 plugin(0),
						 pjets(0) {
  pjets= new vector<fastjet::PseudoJet>();
  vector<fastjet::PseudoJet> particles;
  for( UInt_t i= 0; i < vtl.size(); i++ ) {
    fastjet::PseudoJet pj( vtl[i] );
    if( vindx==0 ) {
      pj.set_user_index( i );
    }
    else {
      pj.set_user_index( vindx->at(i) );
    }
    particles.push_back( pj );
  }
  static map<string,fastjet::JetAlgorithm> jamap;
  if( jamap.size() == 0 ) {
    jamap["kt"]= fastjet::kt_algorithm;
    jamap["cambridge"]= fastjet::cambridge_algorithm;
    jamap["antikt"]= fastjet::antikt_algorithm;
    jamap["genkt"]= fastjet::genkt_algorithm;
    jamap["siscone"]= fastjet::plugin_algorithm;
    jamap["eekt"]= fastjet::ee_kt_algorithm;
    jamap["jade"]= fastjet::plugin_algorithm;
    jamap["eeantikt"]= fastjet::ee_genkt_algorithm;
  }
  string jetalgString( jetalg );
  fastjet::JetAlgorithm ja= jamap[jetalgString];
  //  cout << "TFastJet::TFastJet: " << jetalg << " enum " << ja << endl;
  fastjet::JetDefinition jetdef;
  fastjet::JetDefinition::Recombiner* recombiner= 0;
  if( ja == 99 ) {
    if( jetalgString == "siscone" ) {
      plugin= new fastjet::SISConePlugin( R, 0.75 );
    }
    else if( jetalgString == "jade" ) {
      plugin= new fastjet::JadePlugin();
      recombiner= new EEE0Recombiner();
    }
    else {
      cout << "TFastJet::TFastJet: jet plugin not known: " << ja << endl;
      return;
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
  clusseq= new fastjet::ClusterSequence( particles, jetdef );
}

TFastJet::~TFastJet() {
  if( clusseq ) delete clusseq;
  if( pjets ) delete pjets;
  if( plugin ) delete plugin;
}

const vector<TLorentzVector>& TFastJet::inclusive_jets( const double& ptmin ) {
  vector<fastjet::PseudoJet> incljets= clusseq->inclusive_jets( ptmin );
  *pjets= sorted_by_pt( incljets );
  return copyPseudoJetsToLorentzVectors();
}

const vector<TLorentzVector>& TFastJet::exclusive_jets( const int njets ) {
  vector<fastjet::PseudoJet> excljets= clusseq->exclusive_jets( njets );
  *pjets= sorted_by_E( excljets );
  return copyPseudoJetsToLorentzVectors();
}

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

double TFastJet::ymerge( int njets ) {
  return clusseq->exclusive_ymerge_max( njets );
}

int TFastJet::njets( double ycut ) {
  return clusseq->n_exclusive_jets_ycut( ycut );
}

double TFastJet::Evis() {
  return clusseq->Q();
}

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

