
#include "ObsThrust.hh"
#include "NtupleReader.hh"
#include "DifferentialDataStructure.hh"
#include "MatrixDataStructure.hh"
#include "FilledObservable.hh"
#include "TMath.h"
#include <iostream>
using std::cout;
using std::endl;

ObsThrust::ObsThrust( const vector<Double_t>& bins, 
		      const vector<Analysis>& variations ) :
  ObsDifferential( "thrust", bins ) {
  addAnalyses( variations );
  cout << "ObsThrust::ObsThrust: create " << getName() 
       << " with 1-T and (1-T)^2 weighted distributions" << endl;
  printBinedges();
}

void ObsThrust::addAnalyses( const vector<Analysis>& variations ) {
  // Unweighted thrust and migration matrices
  ObsDifferential::addAnalyses( variations );
  for( size_t ivar= 0; ivar < variations.size(); ivar++ ) {      
    string tag= variations[ivar].getTag();
    weighted1[tag]= new DifferentialDataStructure( binedges );
    weighted2[tag]= new DifferentialDataStructure( binedges );
  }
}

void ObsThrust::fill( NtupleReader* ntr, const Analysis& variation ) {
  string tag= variation.getTag();  
  Double_t thrustvalue= ntr->getThrust( variation.getReco() );
  getAndFillDifferentialDataStructure( thrustvalue, tag, datastructures );
  if( thrustvalue >= 0.0 ) {
    getAndFillDifferentialDataStructure( thrustvalue, tag, weighted1, thrustvalue );
    getAndFillDifferentialDataStructure( thrustvalue, tag, weighted2, 
					 TMath::Power( thrustvalue, 2 ) );
  }
  if( variation.getReco2() != "none" and ntr->isMC() ) {
    Double_t MCthrustvalue= ntr->getThrust( variation.getReco2() );
    matrices[tag]->fill( MCthrustvalue, thrustvalue );
  }
}

vector<FilledObservable*> ObsThrust::getFilledObservables() const {
  cout << "ObsThrust::getFilledObservables: " 
       << "create FilledObservables " << name << " and " << name+"W1" << endl;  
  vector<FilledObservable*> vfobs;
  vfobs.push_back( new FilledObservable( name, datastructures, matrices ) );
  vfobs.push_back( new FilledObservable( name+"W1", weighted1 ) );
  vfobs.push_back( new FilledObservable( name+"W2", weighted2 ) );
  return vfobs;
}

