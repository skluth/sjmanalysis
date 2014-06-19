
#include "ObsThrust.hh"
#include "NtupleReader.hh"
#include "MatrixDataStructure.hh"

ObsThrust::ObsThrust( const vector<Double_t>& bins, 
		      const vector<Analysis>& variations ) :
  ObsDifferential( "thrust", bins ) {
  addAnalyses( variations );
}

void ObsThrust::fill( NtupleReader* ntr, const Analysis& variation ) {
  string tag= variation.getTag();  
  Double_t thrustvalue= ntr->getThrust( variation.getReco() );
  getAndFillDifferentialDataStructure( thrustvalue, tag, datastructures );
  if( variation.getReco2() != "none" and ntr->isMC() ) {
    Double_t MCthrustvalue= ntr->getThrust( variation.getReco2() );
    matrices[tag]->fill( MCthrustvalue, thrustvalue );
  }
}

