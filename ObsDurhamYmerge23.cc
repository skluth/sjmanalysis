
#include "ObsDurhamYmerge23.hh"
#include "NtupleReader.hh"
#include "TMath.h"

ObsDurhamYmerge23::ObsDurhamYmerge23( const vector<Double_t>& bins, 
				      const vector<Analysis>& variations ) :
  ObsDifferential( "durhamymerge23", bins ) {
  addAnalyses( variations );
}

void ObsDurhamYmerge23::fill( NtupleReader* ntr, const Analysis& variation ) {
  Double_t yflip23= ntr->getYmergeD( variation.getReco(), 2 );
  getAndFillDifferentialDataStructure( -TMath::Log10( yflip23 ), variation.getTag(), 
				       datastructures );
}

