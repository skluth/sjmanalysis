
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

// Underflow goes in value[0], overflow goes in values[n+1]:
void DifferentialDataStructure::fill( Double_t value, Double_t weight ) {
  Ntotal++;
  vector<double>::iterator iter= lower_bound( binedges.begin(), binedges.end(), value );
  size_t index= iter-binedges.begin();
  values[index]+= weight;
  errors[index]= TMath::Sqrt( TMath::Power( errors[index], 2 ) + TMath::Power( weight, 2 ) );
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
  Double_t sum= 0.0;
  for( size_t i= 0; i < n-1; i++ ) {
    sum+= values[i+1];
    cout << binedges[i] << " " << binedges[i+1] << " " << values[i+1] << " " 
	 << errors[i+1] << endl;
  }
  cout << "Sum of bins: " << sum << endl;
}

