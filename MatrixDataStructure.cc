
#include "MatrixDataStructure.hh"

using std::vector;

#include <iostream>
using std::cout;
using std::endl;
#include <sstream>
using std::ostringstream;

#include <stdexcept>

// Ctor:
MatrixDataStructure::MatrixDataStructure( const vector<Double_t>& bins ) : 
  binedges(bins), Ntotal(0) {
  ndim= binedges.size()+1;
  array= new Double_t[ndim*ndim];
  for( size_t i= 0; i < ndim; i++ ) {
    for( size_t j= 0; j < ndim; j++ ) {
      array[j*ndim+i]= 0.0;
    }
  }
}

MatrixDataStructure* MatrixDataStructure::clone() const {
  return new MatrixDataStructure( binedges );
}

// Underflows go in value[0], overflow goes in values[n] with n # of binedges:
void MatrixDataStructure::fill( Double_t xvalue, Double_t yvalue ) {
  Ntotal++;
  vector<double>::iterator xiter= lower_bound( binedges.begin(), binedges.end(), 
					       xvalue );
  vector<double>::iterator yiter= lower_bound( binedges.begin(), binedges.end(), 
					       yvalue );
  size_t xindex= xiter-binedges.begin();
  size_t yindex= yiter-binedges.begin();
  array[yindex*ndim+xindex]++;
}

Double_t MatrixDataStructure::getElement( size_t irow, size_t icol ) const {
  checkIndices( irow, icol );
  return array[icol*ndim+irow];
}
void MatrixDataStructure::setElement( size_t irow, size_t icol,
				      Double_t value ) {
  checkIndices( irow, icol );
  array[icol*ndim+irow]= value;
}

void MatrixDataStructure::checkIndices( size_t irow, size_t icol ) const {
  if( irow >= ndim or irow < 0 or icol >= ndim or icol < 0 ) {
    ostringstream txt;
    txt << "Bad irow, icol or ndim: " << irow << " " << icol << " " << ndim;
    throw std::runtime_error( txt.str() );
  }
}

void MatrixDataStructure::Print() const {
  cout << "MatrixDataStructure::Print: " << Ntotal << " events, " 
       << ndim-1 << " binedges:" << endl;
  for( size_t i= 0; i < ndim-1; i++ ) cout << binedges[i] << " ";
  cout << endl;
  for( size_t icol= 0; icol < ndim; icol++ ) {
    for( size_t irow= 0; irow < ndim; irow++ ) {
      cout << array[icol*ndim+irow] << " ";
    }  
    cout << endl;
  }
}

