
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
  Observable( name ), selfCorrelation(scOpt), binedges( bins ),
  nevents(0), neventsb(0), neventsudsc(0) {
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
  datab[tag]= new DifferentialDataStructure( binedges );
  dataudsc[tag]= new DifferentialDataStructure( binedges );
}

// Calculuate EEC as 1/sigma*dEEC/dchi with chi in radian
// incl. self-correlation or not
// weight calculation can be overidden in subclasses
void ObsEEC::fill( NtupleReader* ntr, const Analysis & variation ) {
  const vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( variation.getReco() );
  Int_t iflavour= ntr->getPrimaryFlavour();
  Double_t evis= Evis( vtlv );
  Double_t evis2= evis*evis;
  string tag= variation.getTag();
  DifferentialDataStructure* dds= data.at( tag );
  DifferentialDataStructure* ddsb= datab.at( tag );
  DifferentialDataStructure* ddsudsc= dataudsc.at( tag );
  for( const TLorentzVector& tlv1 : vtlv ) {
    for( const TLorentzVector& tlv2 : vtlv ) {
      if( not selfCorrelation and ( &tlv1 == &tlv2 ) ) continue;
      Double_t angle, eecweight;
      calcWeight( tlv1, tlv2, angle, eecweight );
      Double_t eecvalue= eecweight/evis2;
      // Inclusive
      dds->fill( angle, eecvalue );
      // b or udsc
      if( iflavour == 5 ) {
	ddsb->fill( angle, eecvalue );
      }
      else if( iflavour >= 1 and iflavour <= 4 ) {
	ddsudsc->fill( angle, eecvalue );
      }
    }
  }
  nevents+= 1;
  dds->setNEvents( nevents );
  if( iflavour == 5 ) {
    neventsb+= 1;
    ddsb->setNEvents( neventsb );
  }
  else if( iflavour >= 1 and iflavour <= 4 ) {
    neventsudsc+= 1;
    ddsudsc->setNEvents( neventsudsc );
  }
  return;
}
void ObsEEC::calcWeight( const TLorentzVector& tlv1, const TLorentzVector& tlv2,
			 Double_t& angle, Double_t& weight ) {
  TVector3 v1= tlv1.Vect();
  TVector3 v2= tlv2.Vect();
  angle= v1.Angle( v2 );
  weight= tlv1.E()*tlv2.E();
  return;
}

vector<FilledObservable*> ObsEEC::getFilledObservables() const {
  vector<FilledObservable*> vfobs;
  vfobs.push_back( new FilledObservable( name, data ) );
  vfobs.push_back( new FilledObservable( name+"b", datab ) );
  vfobs.push_back( new FilledObservable( name+"udsc", dataudsc ) );
  return vfobs;
}

