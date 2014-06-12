
#include "JetrateDataStructure.hh"
#include "TMath.h"
#include <iostream>
using std::cout;
using std::endl;

JetrateDataStructure::JetrateDataStructure( const vector<Double_t>& p ) :
  DataStructure() {
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

DataStructure* JetrateDataStructure::clone() {
  return new JetrateDataStructure( points );
}

void JetrateDataStructure::fill( const vector<Double_t>& NJets, Int_t Jetrate ) {
  Ntotal++;
  for( size_t i= 0; i < NJets.size(); i++ ) {
    if( NJets[i] == Jetrate ) {
      values[i]++;
      errors[i]= TMath::Sqrt( TMath::Power( errors[i], 2 ) + 1 );
    }
  }
}

void JetrateDataStructure::normalise() {
  for( size_t i= 0; i < values.size(); i++ ) {
    values[i]/= Ntotal;
    errors[i]= TMath::Sqrt( values[i]*(1.0-values[i])/Ntotal ); 
  }
}

void JetrateDataStructure::print() {
  size_t n= points.size();
  for( size_t i= 0; i < n; i++ ) {
    cout << points[i] << " " << values[i] << " " << errors[i] << endl;
  }
}

