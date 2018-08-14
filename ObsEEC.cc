
#include "ObsEEC.hh"

#include "NtupleReader.hh"
#include "DifferentialDataStructure.hh"
#include "FilledObservable.hh"

#include "TLorentzVector.h"
#include "TVector3.h"

#include <iostream>
using std::cout;
using std::endl;

using std::string;
using std::vector;

ObsEEC::ObsEEC( const string & name,
		const vector<Double_t> & bins,
		const vector<Analysis> & variations,
		const bool scOpt,
		const bool lprint ) :
  Observable( name ), binedges( bins ), selfCorrelation(scOpt) {
  addAnalyses( variations );
  if( lprint ) {
    cout << "ObsEEC::ObsEEC: Energy-Energy Correlation for " << name;
    if( not selfCorrelation ) cout << " w/o self-correlation";
    cout << endl;
    printVectorD( "Binedges:", bins );
  }
  return;
}

ObsEEC::~ObsEEC() {}

void ObsEEC::addAnalysis( const Analysis & analysis ) {
  string tag= analysis.getTag();
  data[tag]= new DifferentialDataStructure( binedges );
}

// Calculuate EEC as 1/sigma*dEEC/dchi with chi in radian
// incl. self-correlation or not
void ObsEEC::fill( NtupleReader* ntr, const Analysis & variation ) {
  const vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( variation.getReco() );
  Double_t evis= ntr->Evis( vtlv );
  Double_t evis2= evis*evis;
  string tag= variation.getTag();
  DifferentialDataStructure* dds= data.at( tag );
  Double_t nevents= dds->getNEvents();
  for( const TLorentzVector& tlv1 : vtlv ) {
    for( const TLorentzVector& tlv2 : vtlv ) {
      if( not selfCorrelation and ( &tlv1 == &tlv2 ) ) continue;
      TVector3 v1= tlv1.Vect();
      TVector3 v2= tlv2.Vect();
      Double_t angle= v1.Angle( v2 );
      Double_t eecvalue= tlv1.E()*tlv2.E()/evis2;
      dds->fill( angle, eecvalue );
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

