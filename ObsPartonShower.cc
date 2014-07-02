
#include "ObsPartonShower.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
#include "DifferentialDataStructure.hh"
#include "MatrixDataStructure.hh"
#include "FilledObservable.hh"
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>
using std::cout;
using std::endl;

ObsPartonShower::ObsPartonShower( const vector<Double_t>& a14bins,
				  const vector<Double_t>& c202bins,
				  const vector<Double_t>& asbins,
				  const vector<Analysis>& variations,
				  Double_t y34c ) :
  ObsDifferential( "partonshower", a14bins ), 
  c202binedges(c202bins), asbinedges(asbins), y34cut(y34c) {
  addAnalyses( variations );
  cout << "ObsPartonShower::ObsPartonShower: create " << getName() 
       << " with A14, C202, AS with cut y34cut=" << y34cut << endl;
  printVectorD( "A14 binedges:", a14bins );
  printVectorD( "C202 binedges:", c202bins );
  printVectorD( "AS binedges:", asbins );
}

void ObsPartonShower::addAnalyses( const vector<Analysis>& variations ) {
  // A14
  ObsDifferential::addAnalyses( variations );
  // C202 and AS
  for( size_t ivar= 0; ivar < variations.size(); ivar++ ) {      
    string tag= variations[ivar].getTag();
    c202data[tag]= new DifferentialDataStructure( c202binedges );
    asdata[tag]= new DifferentialDataStructure( asbinedges );
  }
}


void ObsPartonShower::fill( NtupleReader* ntr, const Analysis& variation ) {
  string reco= variation.getReco();
  //  Double_t a14value= getA14value( ntr, reco );
  vector<Double_t> values= getValues( ntr, reco );
  string tag= variation.getTag();
  getAndFillDifferentialDataStructure( values[0], tag, datastructures );
  getAndFillDifferentialDataStructure( values[1], tag, c202data );
  getAndFillDifferentialDataStructure( values[2], tag, asdata );
  string reco2= variation.getReco2();
  if( reco2 != "none" and ntr->isMC() ) {
    //    Double_t MCa14value= getA14value( ntr, reco2 );
    vector<Double_t> MCvalues= getValues( ntr, reco2 );    
    matrices[tag]->fill( MCvalues[0], values[0] );
  }
}

vector<Double_t> ObsPartonShower::getValues( NtupleReader* ntr, 
						    const string& reco ) {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( reco );
  TFastJet tfj( vtlv, "eekt" );
  //  Double_t a14value= -1.0;
  vector<Double_t> values( 3 );
  for( size_t i= 0; i < 3; i++ ) values[i]= -1.0;
  if( tfj.ymerge( 3 ) > y34cut ) {
    vector<TLorentzVector> jets= tfj.exclusive_jets( 4 );
    Double_t angle12= calcAngle( jets[0], jets[1] );
    Double_t angle13= calcAngle( jets[0], jets[2] );
    Double_t angle23= calcAngle( jets[1], jets[2] );
    if( angle12 > 2.0*TMath::Pi()/3.0 and angle13 > 2.0*TMath::Pi()/3.0 and
	angle23 < TMath::Pi()/6.0 ) { 
      values[0]= calcAngle( jets[0], jets[3] )/TMath::Pi(); // A14
      Double_t angle24= calcAngle( jets[1], jets[3] );
      if( angle24 < TMath::Pi()/2.0 ) {
	values[2]= (angle24-angle23)/TMath::Pi(); // AS
      }
    }
  }
  return values;
}

Double_t ObsPartonShower::calcAngle( const TLorentzVector& tlv1, 
				     const TLorentzVector& tlv2 ) {
  const TVector3& tv1= tlv1.Vect();
  const TVector3& tv2= tlv2.Vect();
  Double_t angle= tv1.Angle( tv2 );
  return angle;
}

vector<FilledObservable*> ObsPartonShower::getFilledObservables() const {
  cout << "ObsPartonShower::getFilledObservables: " 
       << "create FilledObservables a14, c202, as" << endl;  
  vector<FilledObservable*> vfobs;
  vfobs.push_back( new FilledObservable( "a14", datastructures, matrices ) );
  vfobs.push_back( new FilledObservable( "c202", c202data ) );
  vfobs.push_back( new FilledObservable( "as", asdata ) );
  return vfobs;
}
