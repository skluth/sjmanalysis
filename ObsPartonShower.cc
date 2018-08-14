
#include "ObsPartonShower.hh"

#include "NtupleReader.hh"
#include "TFastJet.hh"
#include "DifferentialDataStructure.hh"
#include "MatrixDataStructure.hh"
#include "FilledObservable.hh"
#include "TLorentzVector.h"

#include "TMath.h"

#include <sstream>
using std::ostringstream;
#include <iostream>
using std::cout;
using std::endl;
#include <stdexcept>
using std::logic_error;

using std::string;
using std::vector;

ObsPartonShower::ObsPartonShower( const vector<Double_t>& a14bins,
				  const vector<Double_t>& c202bins,
				  const vector<Double_t>& asbins,
				  const vector<Double_t>& mrbins,
				  const vector<Analysis>& variations,
				  Double_t y34c, Double_t y34y23c,
				  bool lprint ) :
  Observable( "partonshower" ),
  a14binedges(a14bins), c202binedges(c202bins), asbinedges(asbins), mrbinedges(mrbins),
  y34cut(y34c), y34y23cut(y34y23c) {
  addAnalyses( variations );
  if( lprint ) {
    cout << "ObsPartonShower::ObsPartonShower: create " << getName() 
	 << " with A14, C202, AS, MR with cuts y34cut=" << y34cut 
	 << " and y34y23cut=" << y34y23cut << endl;
    printVectorD( "A14 binedges:", a14bins );
    printVectorD( "C202 binedges:", c202bins );
    printVectorD( "AS binedges:", asbins );
    printVectorD( "MR binedges:", mrbins );
  }
}

//void ObsPartonShower::addAnalyses( const vector<Analysis>& variations ) {
void ObsPartonShower::addAnalysis( const Analysis& analysis ) {
  // A14, C202, AS, MR:
  string tag= analysis.getTag();
  a14data[tag]= new DifferentialDataStructure( a14binedges );
  c202data[tag]= new DifferentialDataStructure( c202binedges );
  asdata[tag]= new DifferentialDataStructure( asbinedges );
  mrdata[tag]= new DifferentialDataStructure( mrbinedges );
  if( analysis.getReco2() != "none" ) {
    a14matrices[tag]= new MatrixDataStructure( a14binedges );
    c202matrices[tag]= new MatrixDataStructure( c202binedges );
    asmatrices[tag]= new MatrixDataStructure( asbinedges );
    mrmatrices[tag]= new MatrixDataStructure( mrbinedges );
  }
  return;
}

void ObsPartonShower::fill( NtupleReader* ntr, const Analysis& variation ) {
  string reco= variation.getReco();
  string tag= variation.getTag();
  vector<Double_t> values= getValues( ntr, reco );
  a14data.at(tag)->fill( values[0] );
  c202data.at(tag)->fill( values[1] );
  asdata.at(tag)->fill( values[2] );
  mrdata.at(tag)->fill( values[3] );
  string reco2= variation.getReco2();
  if( reco2 != "none" and ntr->isMC() ) {
    vector<Double_t> MCvalues= getValues( ntr, reco2 );    
    a14matrices.at(tag)->fill( MCvalues[0], values[0] );
    c202matrices.at(tag)->fill( MCvalues[1], values[1] );
    asmatrices.at(tag)->fill( MCvalues[2], values[2] );
    mrmatrices.at(tag)->fill( MCvalues[3], values[3] );
  }
  return;
}

// From Nadines code example A14.C, C202.C, AS.C, MR.C
vector<Double_t> ObsPartonShower::getValues( NtupleReader* ntr, 
					     const string& reco ) {
  const vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( reco );
  TFastJet tfj( vtlv, "eekt" );
  vector<Double_t> values{ -1.0, -1.0, -1.0, -1.0 };
  if( tfj.ymerge( 3 ) > y34cut ) {
    vector<TLorentzVector> jets= tfj.exclusive_eejets( 4 );
    Double_t angle12= calcAngle( jets[0], jets[1] );
    Double_t angle13= calcAngle( jets[0], jets[2] );
    Double_t angle23= calcAngle( jets[1], jets[2] );
    if( angle12 > 2.0*TMath::Pi()/3.0 and angle13 > 2.0*TMath::Pi()/3.0 and
	angle23 < TMath::Pi()/6.0 ) { 
      values[0]= calcAngle( jets[0], jets[3] )/TMath::Pi(); // A14
      values[1]= 
	EnergyCorrelator( jets, 1, 0.2 ) * 
	EnergyCorrelator( jets, 3, 0.2 ) /
	TMath::Power( EnergyCorrelator( jets, 2, 0.2 ), 2 ); // C202
      Double_t angle24= calcAngle( jets[1], jets[3] );
      if( angle24 < TMath::Pi()/2.0 ) {
	values[2]= ( angle24 - angle23 )/TMath::Pi(); // AS
      }
    }
    if( tfj.ymerge( 3 )/tfj.ymerge( 2 ) > y34y23cut ) {
      jets= tfj.exclusive_eejets( 2 );
      values[3]= TMath::Min( jets[1].M2(), jets[0].M2() ) /
	TMath::Max( jets[1].M2(), jets[0].M2() ); // MR
    }
  }
  return values;
}

Double_t ObsPartonShower::calcAngle( const TLorentzVector& tlv1, 
				     const TLorentzVector& tlv2 ) {
  TVector3 tv1= tlv1.Vect();
  TVector3 tv2= tlv2.Vect();
  Double_t angle= tv1.Angle( tv2 );
  return angle;
}

Double_t ObsPartonShower::EnergyCorrelator( const vector<TLorentzVector>& jets, 
					    Int_t N, Double_t beta ) {
  Double_t answer= 0.0;
  if( N == 1 ) {
    for( size_t i= 0; i < jets.size(); i++ ) {
      answer+= jets[i].E();
    }
  }
  else if( N == 2 ) {
    for( size_t i= 0; i < jets.size(); i++ ) {
      for( size_t j= i+1; j < jets.size(); j++ ) { 
	answer+= jets[i].E() * jets[j].E()
	  * TMath::Power( calcAngle( jets[i], jets[j] ), beta );
      }
    }
  }
  else if( N == 3 ) {
    for( size_t i= 0; i < jets.size(); i++ ) {
      for( size_t j= i+1; j < jets.size(); j++ ) {
	Double_t ans_ij= jets[i].E() * jets[j].E()
	  * TMath::Power( calcAngle( jets[i], jets[j] ), beta );
	for( size_t k= j+1; k < jets.size(); k++ ) {
	  answer+= ans_ij * jets[k].E()
	    * TMath::Power( calcAngle( jets[j], jets[k] ), beta )
	    * TMath::Power( calcAngle( jets[i], jets[k] ), beta );
	}
      }
    } 
  }
  else {
    ostringstream txt;
    txt << "ObsPartonShower::EnergyCorrelator: wrong N: " << N;
    throw logic_error( txt.str() );
  }
  return answer;
}

vector<FilledObservable*> ObsPartonShower::getFilledObservables() const {
  // cout << "ObsPartonShower::getFilledObservables: " 
  //      << "create FilledObservables a14, c202, as, mr" << endl;  
  vector<FilledObservable*> vfobs;
  vfobs.push_back( new FilledObservable( "a14", a14data, a14matrices ) );
  vfobs.push_back( new FilledObservable( "c202", c202data, c202matrices ) );
  vfobs.push_back( new FilledObservable( "as", asdata, asmatrices ) );
  vfobs.push_back( new FilledObservable( "mr", mrdata, mrmatrices ) );
  return vfobs;
}


