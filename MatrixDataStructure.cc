
#include "MatrixDataStructure.hh"

#include <algorithm>
using std::upper_bound;

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
      array[i*ndim+j]= 0.0;
    }
  }
}

MatrixDataStructure* MatrixDataStructure::clone() const {
  return new MatrixDataStructure( binedges );
}

// Underflows go in value[0], overflow goes in values[n] with n # of binedges:
void MatrixDataStructure::fill( Double_t xvalue, Double_t yvalue ) {
  Ntotal++;
  vector<double>::iterator xiter= upper_bound( binedges.begin(), binedges.end(), 
					       xvalue );
  vector<double>::iterator yiter= upper_bound( binedges.begin(), binedges.end(), 
					       yvalue );
  size_t icol= xiter-binedges.begin();
  size_t irow= yiter-binedges.begin();
  array[icol*ndim+irow]++;
}

Double_t MatrixDataStructure::getElement( size_t icol, size_t irow ) const {
  checkIndices( icol, irow );
  return array[icol*ndim+irow];
}
void MatrixDataStructure::setElement( size_t icol, size_t irow,
				      Double_t value ) {
  checkIndices( icol, irow );
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

// Normalise columns w/o under- and overflow
void MatrixDataStructure::normaliseColumns() {
  for( size_t icol= 1; icol < ndim-1; icol++ ) {
    Double_t sumcol= 0.0;
    for( size_t irow= 1; irow < ndim-1; irow++ ) sumcol+= array[icol*ndim+irow];
    if( sumcol != 0.0 ) {
      for( size_t irow= 1; irow < ndim-1; irow++ ) array[icol*ndim+irow]/= sumcol;
    }    
  }
}

MatrixDataStructure* applyEfficiency( const MatrixDataStructure & m,
				      const vector<Double_t> & v ) {
  size_t mndim= m.getNdim();
  size_t vndim= v.size();
  if( mndim != vndim ) {
    throw std::logic_error( "Matrix and vector dimensions don't match" );
  }
  MatrixDataStructure* result= new MatrixDataStructure( m.getBinedges() );
  for( size_t icol= 0; icol < mndim; icol++ ) {
    for( size_t irow= 0; irow < mndim; irow++ ) {
      result->setElement( icol, irow, m.getElement( icol, irow )*v[irow] );
    }
  }
  return result;
}

vector<Double_t> multiply( const MatrixDataStructure & m,
			   const vector<Double_t> & v ) {
  size_t mndim= m.getNdim();
  size_t vndim= v.size();
  if( mndim != vndim ) {
    throw std::logic_error( "Matrix and vector dimensions don't match" );
  }
  vector<Double_t> result( vndim );
  for( size_t irow= 0; irow < mndim; irow++ ) {
    Double_t sumrow= 0.0;
    for( size_t icol= 0; icol < mndim; icol++ ) {
      sumrow+= m.getElement( icol, irow )*v[icol];
    }
    result[irow]= sumrow;
  }
  return result;
}

MatrixDataStructure* similarityVector( const MatrixDataStructure & m,
				       const vector<Double_t> & v ) {
  size_t mndim= m.getNdim();
  size_t vndim= v.size();
  if( mndim != vndim ) {
    throw std::logic_error( "Matrix and vector dimensions don't match" );
  }
  MatrixDataStructure* result= new MatrixDataStructure( m.getBinedges() );
  for( size_t icol= 0; icol < mndim; icol++ ) {
    for( size_t irow= 0; irow < mndim; irow++ ) {
      Double_t sumk= 0.0;
      for( size_t k= 0; k < mndim; k++ ) {
	sumk+= v[k]*v[k]*m.getElement( k, icol )*m.getElement( k, irow );
      }
      result->setElement( icol, irow, sumk );
    }
  }
  return result;
}
