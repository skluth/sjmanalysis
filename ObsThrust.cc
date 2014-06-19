
#include "ObsThrust.hh"
#include "NtupleReader.hh"
#include "MatrixDataStructure.hh"

#include <iostream>
using std::cout;
using std::endl;
#include "TMath.h"

ObsThrust::ObsThrust( const vector<Double_t>& bins, 
		      const vector<Analysis>& variations ) :
  ObsDifferential( "thrust", bins ) {
  addAnalyses( variations );

  for( size_t ivar= 0; ivar < variations.size(); ivar++ ) {      
    if( variations[ivar].getReco2() != "none" ) {
      string tag= variations[ivar].getTag();
      matrices[tag]= new MatrixDataStructure( binedges );
    }
  }

}

void ObsThrust::fill( NtupleReader* ntr, const Analysis& variation ) {
  string tag= variation.getTag();  
  Double_t thrustvalue= ntr->getThrust( variation.getReco() );
  getAndFillDifferentialDataStructure( thrustvalue, tag, datastructures );
  if( variation.getReco2() != "none" and ntr->isMC() ) {
    Double_t MCthrustvalue= ntr->getThrust( variation.getReco2() );
    MatrixDataStructure* matrix= matrices[tag];
    matrix->fill( MCthrustvalue, thrustvalue );
  }
}

