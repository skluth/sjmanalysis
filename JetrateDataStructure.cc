
#include "JetrateDataStructure.hh"
#include "MatrixDataStructure.hh"

#include "TMath.h"
#include <iostream>
using std::cout;
using std::endl;
#include <sstream>
using std::ostringstream;
#include <stdexcept>
using std::logic_error;

using std::vector;

JetrateDataStructure::JetrateDataStructure( const vector<Double_t>& p, 
					    Int_t njet ) :
  DataStructure(), Jetrate(njet) {
  size_t n= p.size();
  points.resize( n );
  values.resize( n );
  errors.resize( n );
  for( size_t i= 0; i < n; i++ ) {
    points[i]= p[i];
    values[i]= 0.0;
    errors[i]= 0.0;
  }
}

DataStructure* JetrateDataStructure::clone() const {
  return new JetrateDataStructure( points, Jetrate );
}

// Set error matrix for tag with bins for point numbers:
void JetrateDataStructure::setErrorMatrix() {
  vector<Double_t> bins;
  for( size_t i= 0; i <= points.size(); i++ ) bins.push_back( i+0.5 );
  errorMatrix= new MatrixDataStructure( bins );
}

void JetrateDataStructure::fill( const vector<Double_t>& NJets ) {
  if( NJets.size() != values.size() ) {
    ostringstream txt;
    txt << "Input vector wrong size " << NJets.size() 
	<< ", expected " << values.size();
    throw logic_error( txt.str() );
  }
  Ntotal++;
  for( size_t i= 0; i < NJets.size(); i++ ) {
    if( NJets[i] == Jetrate ) {
      values[i]++;
      errors[i]= TMath::Sqrt( TMath::Power( errors[i], 2 ) + 1.0 );
    }
  }
}

void JetrateDataStructure::normalise() {
  checkNormalised();
  checkNtotalGTZero();
  for( size_t i= 0; i < values.size(); i++ ) {
    values[i]/= Ntotal;
    errors[i]= TMath::Sqrt( values[i]*(1.0-values[i])/Ntotal ); 
  }
  setNormalisedTrue();
}

void JetrateDataStructure::Print() const {
  size_t n= points.size();
  for( size_t i= 0; i < n; i++ ) {
    cout << points[i] << " " << values[i] << " " << errors[i] << endl;
  }
}

