
#include "MatrixDataStructure.hh"
#include <iostream>
using std::cout;
using std::endl;

MatrixDataStructure::MatrixDataStructure( vector<Double_t> bins ) : 
  binedges(bins), Ntotal(0) {
  ndim= binedges.size()+1;
  array= new Double_t[ndim*ndim];
}

// Underflows go in value[0], overflow goes in values[n] with n # of binedges:
void MatrixDataStructure::fill( Double_t xvalue, Double_t yvalue ) {
  Ntotal++;
  vector<double>::iterator xiter= lower_bound( binedges.begin(), binedges.end(), xvalue );
  vector<double>::iterator yiter= lower_bound( binedges.begin(), binedges.end(), yvalue );
  size_t xindex= xiter-binedges.begin();
  size_t yindex= yiter-binedges.begin();
  array[yindex*ndim+xindex]++;
}

Double_t MatrixDataStructure::getElement( size_t irow, size_t icol ) const {
  if( irow >= ndim ) cout << "irow > ndim " << irow << " " << ndim << endl;
  if( irow < 0 ) cout << "irow < 0 " << irow << endl;
  if( icol >= ndim ) cout << "icol > ndim " << icol << " " << ndim << endl;
  if( icol < 0 ) cout << "icol < 0 " << icol << endl;
  return array[icol*ndim+irow];
}

void MatrixDataStructure::print() const {
  cout << "MatrixDataStructure::print: " << Ntotal << " events, " 
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

