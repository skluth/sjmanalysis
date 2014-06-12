
#include "DifferentialDataStructure.hh"
#include "TMath.h"
#include <algorithm>
using std::lower_bound;
#include <iostream>
using std::cout;
using std::endl;

DifferentialDataStructure::DifferentialDataStructure( const vector<Double_t>& bins ) : 
  DataStructure(), binedges(bins) {
  size_t n= binedges.size();
  points.resize( n-1 );
  values.resize( n+1 );
  errors.resize( n+1 );
  for( size_t i= 0; i < n-1; i++ ) {
    points[i]= (bins[i+1]+bins[i])/2.0;
  }
  for( size_t i= 0; i < n+1; i++ ) {
    values[i]= 0.0;
    errors[i]= 0.0;
  }
}

DataStructure* DifferentialDataStructure::clone() {
  return new DifferentialDataStructure( binedges );
}

// Return w/o under- and overflows, i.e. proper data points only:
vector<Double_t> DifferentialDataStructure::getValues() const { 
  vector<Double_t> result( values.begin()+1, values.end()-1 );
  return result; 
}
vector<Double_t> DifferentialDataStructure::getErrors() const { 
  vector<Double_t> result( errors.begin()+1, errors.end()-1 );
  return result;
}

void DifferentialDataStructure::setValues( const vector<Double_t>& valuesin ) {
  values[0]= 0.0;
  size_t n= valuesin.size();
  values[n+2]= 0.0;
  for( size_t i= 0; i < n; i++ ) {
    values[i+1]= valuesin[i];
  }    
} 

void DifferentialDataStructure::setErrors( const vector<Double_t>& errorsin ) {
  errors[0]= 0.0;
  size_t n= errorsin.size();
  errors[n+2]= 0.0;
  for( size_t i= 0; i < n; i++ ) {
    errors[i+1]= errorsin[i];
  }    
} 

// Underflow goes in value[0], overflow goes in values[n+1]:
void DifferentialDataStructure::fill( Double_t value ) {
  Ntotal++;
  vector<double>::iterator iter= lower_bound( binedges.begin(), binedges.end(), value );
  size_t index= iter-binedges.begin();
  values[index]++;
  errors[index]= TMath::Sqrt( TMath::Power( errors[index], 2 ) + 1 );
}

// Normalise only bins, not under- or underflow, because binwidth is not defined:
void DifferentialDataStructure::normalise() {
  if( Ntotal > 0.0 ) {
    for( size_t i= 0; i < points.size(); i++ ) {
      Double_t binw= binedges[i+1]-binedges[i];
      values[i+1]/= (Ntotal*binw);
      errors[i+1]/= (Ntotal*binw);
    }
  }
}

void DifferentialDataStructure::print() {
  size_t n= binedges.size();
  cout << "Under- and overflow: " << values[0] << " +/- " << errors[0] 
       << ", " << values[n] << " +/- " << errors[n] << endl;
  for( size_t i= 0; i < n-1; i++ ) {
    cout << binedges[i] << " " << binedges[i+1] << " " << values[i+1] << " " 
	 << errors[i+1] << endl;
  }
}

