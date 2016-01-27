
#include "ObsDifferential.hh"
#include "NtupleReader.hh"
#include "DataStructure.hh"
#include "DifferentialDataStructure.hh"
#include "MatrixDataStructure.hh"
#include "FilledObservable.hh"
#include "DifferentialCalculator.hh"

#include <iostream>
using std::cout;
using std::endl;

ObsDifferential::ObsDifferential( const string& name,
				  const vector<Double_t>& bins,
				  const vector<Analysis>& variations,
				  const DifferentialCalculator* calc,
				  const bool lprint ) :
  Observable(name), binedges(bins), calculator(calc) {
  addAnalyses( variations );
  // for( size_t ivar= 0; ivar < variations.size(); ivar++ ) {      
  //   string tag= variations[ivar].getTag();
  //   data[tag]= new DifferentialDataStructure( bins );
  //   weighted1[tag]= new DifferentialDataStructure( bins );
  //   weighted2[tag]= new DifferentialDataStructure( bins );
  //   if( variations[ivar].getReco2() != "none" ) {
  //     matrices[tag]= new MatrixDataStructure( bins );
  //   }
  // }
  if( lprint ) {
    cout << "ObsDifferential::ObsDifferential: ds/dy, ds/dy*y, ds/dy*y**2 for " 
	 << name << endl;
    printVectorD( "Binedges:", bins );

  }
  return;
}

ObsDifferential::~ObsDifferential() {}

void ObsDifferential::addAnalyses( const vector<Analysis>& variations ) {
  for( size_t ivar= 0; ivar < variations.size(); ivar++ ) {      
    string tag= variations[ivar].getTag();
    data[tag]= new DifferentialDataStructure( binedges );
    weighted1[tag]= new DifferentialDataStructure( binedges );
    weighted2[tag]= new DifferentialDataStructure( binedges );
    if( variations[ivar].getReco2() != "none" ) {
      matrices[tag]= new MatrixDataStructure( binedges );
    }
  }
}


bool ObsDifferential::containsAnalysis( const Analysis& analysis ) {
  return data.find( analysis.getTag() ) != data.end();
}

void ObsDifferential::fill( NtupleReader* ntr, const Analysis& variation ) {
  string tag= variation.getTag();  
  Double_t value= calculator->getValue( ntr, variation.getReco() );
  data[tag]->fill( value );
  if( value >= 0.0 ) {
    weighted1[tag]->fill( value, value );
    weighted2[tag]->fill( value, TMath::Power( value, 2 ) );
  }
  if( variation.getReco2() != "none" and ntr->isMC() ) {
    Double_t MCvalue= calculator->getValue( ntr, variation.getReco2() );
    matrices[tag]->fill( MCvalue, value );
  }
}

vector<FilledObservable*> ObsDifferential::getFilledObservables() const {
  // cout << "ObsDifferential::getFilledObservables: " 
  //      << "create FilledObservables " << name << " and " << name+"W1" << endl;  
  vector<FilledObservable*> vfobs;
  vfobs.push_back( new FilledObservable( name, data, matrices ) );
  vfobs.push_back( new FilledObservable( name+"W1", weighted1 ) );
  vfobs.push_back( new FilledObservable( name+"W2", weighted2 ) );
  return vfobs;
}

