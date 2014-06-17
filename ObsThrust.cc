
#include "ObsThrust.hh"
#include "NtupleReader.hh"

#include <iostream>
using std::cout;
using std::endl;
#include "TMath.h"

ObsThrust::ObsThrust( const vector<Double_t>& bins, 
		      const vector<Analysis>& variations ) :
  ObsDifferential( "thrust", bins ) {
  addAnalyses( variations );
}

void ObsThrust::fill( NtupleReader* ntr, const Analysis& variation ) {
  Double_t thrustvalue= ntr->getThrust( variation.getReco() );
  getAndFillDifferentialDataStructure( thrustvalue, variation.getTag(), datastructures );
}

