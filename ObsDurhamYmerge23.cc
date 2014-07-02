
#include "ObsDurhamYmerge23.hh"
#include "NtupleReader.hh"
#include "MatrixDataStructure.hh"
#include "TMath.h"
#include <iostream>
using std::cout;
using std::endl;

ObsDurhamYmerge23::ObsDurhamYmerge23( const vector<Double_t>& bins, 
				      const vector<Analysis>& variations ) :
  ObsDifferential( "durhamymerge23", bins ) {
  addAnalyses( variations );
  cout << "ObsDurhamYmerge23::ObsDurhamYmerge23: create " << getName() << endl;
  printBinedges();
}

void ObsDurhamYmerge23::fill( NtupleReader* ntr, const Analysis& variation ) {
  string tag= variation.getTag();
  Double_t log10yflip23= -TMath::Log10( ntr->getYmergeD( variation.getReco(), 2 ) );
  getAndFillDifferentialDataStructure( log10yflip23, tag, datastructures );
  if( variation.getReco2() != "none" and ntr->isMC() ) {
    Double_t MClog10yflip23= -TMath::Log10( ntr->getYmergeD( variation.getReco2(), 2 ) );
    matrices[tag]->fill( MClog10yflip23, log10yflip23 );
  }
}

