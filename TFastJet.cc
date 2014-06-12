
#include "TFastJet.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/SISConePlugin.hh"
#include "fastjet/JadePlugin.hh"

#include "TParticle.h"
#include "TLorentzVector.h"

#include <map>
#include <string>
#include <iostream>

using std::cout;
using std::endl;
using std::string;
using std::map;

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
  if( ja == 99 ) {
    if( jetalgString == "siscone" ) {
      plugin= new fastjet::SISConePlugin( R, 0.75 );
    }
    else if( jetalgString == "jade" ) {
      plugin= new fastjet::JadePlugin();
    }
    else {
      cout << "TFastJet::TFastJet: jet plugin not known: " << ja << endl;
      return;
    }
    jetdef= fastjet::JetDefinition( plugin );
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

vector<TLorentzVector>& TFastJet::inclusive_jets( const double& ptmin ) {
  vector<fastjet::PseudoJet> incljets= clusseq->inclusive_jets( ptmin );
  *pjets= sorted_by_pt( incljets );
  return copyPseudoJetsToLorentzVectors();
}

vector<TLorentzVector>& TFastJet::exclusive_jets( const int njets ) {
  vector<fastjet::PseudoJet> excljets= clusseq->exclusive_jets( njets );
  *pjets= sorted_by_E( excljets );
  return copyPseudoJetsToLorentzVectors();
}

vector<TLorentzVector>& TFastJet::copyPseudoJetsToLorentzVectors() {
  vector<TLorentzVector>* jetstlv= new vector<TLorentzVector>();
  fastjet::PseudoJet pj;
  for( UInt_t i= 0; i < pjets->size(); i++ ) {
    pj= (*pjets)[i];
    TLorentzVector tlv( pj.px(), pj.py(), pj.pz(), pj.E() );
    jetstlv->push_back( tlv );
  }
  return *jetstlv;
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

// g++ -c TFastJet.cc -I ../fastjet/fastjet-2.4.2/include/ -I /usr/include/root/

// g++ -shared -fPIC -o TFastJet.so -I /usr/include/root/ -I ../fastjet/fastjet-2.4.2/install/include/ -L ../fastjet/fastjet-2.4.2/install/lib -lfastjet -lSISConePlugin -lsiscone TFastJet.cc TFastJetDict.cc

// and don't forget to set LD_LIBRARY_PATH for python

