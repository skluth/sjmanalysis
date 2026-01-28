
#include "ObsJEEC.hh"

//#include "NtupleReader.hh"
//#include "DifferentialDataStructure.hh"
//#include "FilledObservable.hh"

#include "TLorentzVector.h"
#include "TVector3.h"

#include <iostream>
using std::cout;
using std::endl;

using std::string;
using std::vector;

#include <cmath>

ObsJEEC::ObsJEEC( const string & name,
		  const vector<Double_t> & bins,
		  const vector<Analysis> & variations,
		  const bool scOpt,
		  const bool lprint ) :
  ObsEEC( name, bins, variations, scOpt, lprint ) {
  addAnalyses( variations );
  if( lprint ) {
    cout << "ObsJEEC::ObsJEEC: JADE Energy-Energy Correlation for " << name;
    if( not selfCorrelation ) cout << " w/o self-correlation";
    cout << endl;
    printVectorD( "Binedges:", bins );
  }
  return;
}

ObsJEEC::~ObsJEEC() {}

void ObsJEEC::calcWeight( const TLorentzVector& tlv1, const TLorentzVector& tlv2,
			  Double_t& angle, Double_t& weight ) {
  TVector3 v1= tlv1.Vect();
  TVector3 v2= tlv2.Vect();
  angle= v1.Angle( v2 );
  weight= tlv1.E()*tlv2.E()*(1.0-cos(angle));
  return;
}

