
#include "ObsEEC.hh"
#include "NtupleReader.hh"
#include "DataStructure.hh"
#include "DifferentialDataStructure.hh"
#include "FilledObservable.hh"
#include "TLorentzVector.h"
#include "TVector3.h"

#include <iostream>
using std::cout;
using std::endl;

ObsEEC::ObsEEC( const vector<Double_t>& bins,
		const vector<Analysis>& variations,
		const bool lprint ) :
  Observable( "EEC" ), binedges( bins ), rad2grad(180.0/3.141592653589793) {
  addAnalyses( variations );
  if( lprint ) {
    cout << "ObsEEC::ObsEEC: Energy-Energy Correlation for " 
	 << name << endl;
    printVectorD( "Binedges:", bins );
  }
  return;
}

ObsEEC::~ObsEEC() {}

void ObsEEC::addAnalysis( const Analysis& analysis ) {
  string tag= analysis.getTag();
  data[tag]= new DifferentialDataStructure( binedges );
}

void ObsEEC::fill( NtupleReader* ntr, const Analysis& variation ) {
  string tag= variation.getTag();  
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( variation.getReco() );
  Double_t evis= ntr->Evis( vtlv );
  Double_t evis2= evis*evis;
  DifferentialDataStructure* dds= data.at( tag );
  Double_t nevents= dds->getNEvents();
  for( const TLorentzVector& tlv1 : vtlv ) {
    for( const TLorentzVector& tlv2 : vtlv ) {
      if( &tlv1 == &tlv2 ) {
	dds->fill( 0.0, tlv1.E()/evis );
      }
      else {
	TVector3 v1= tlv1.Vect();
	TVector3 v2= tlv2.Vect();
	//      Double_t angle= v1.Angle( v2 )*rad2grad;
	Double_t angle= v1.Angle( v2 );
	Double_t eecvalue= tlv1.E()*tlv2.E()/evis2;
	dds->fill( angle, eecvalue );
      }
    }
  }
  nevents+= 1.0;
  dds->setNEvents( nevents );
  return;
}

vector<FilledObservable*> ObsEEC::getFilledObservables() const {
  vector<FilledObservable*> vfobs;
  vfobs.push_back( new FilledObservable( name, data ) );
  return vfobs;
}

