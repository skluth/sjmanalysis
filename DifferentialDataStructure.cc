
#include "DifferentialDataStructure.hh"
#include "MatrixDataStructure.hh"

#include "TMath.h"
#include <algorithm>
using std::upper_bound;
#include <iostream>
using std::cout;
using std::endl;
#include <sstream>
using std::ostringstream;
#include <stdexcept>
using std::logic_error;

using std::vector;

DifferentialDataStructure::DifferentialDataStructure( const vector<Double_t> & bins ) : 
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

DataStructure* DifferentialDataStructure::clone() const {
  return new DifferentialDataStructure( binedges );
}

// Underflow goes in values[0], overflow goes in values[n+1]:
void DifferentialDataStructure::fill( Double_t value, Double_t weight ) {
  Ntotal++;
  vector<double>::iterator iter= upper_bound( binedges.begin(),
					      binedges.end(),
					      value );
  size_t index= iter-binedges.begin();
  values[index]+= weight;
  errors[index]= TMath::Sqrt( TMath::Power( errors[index], 2 ) +
			      TMath::Power( weight, 2 ) );
}

// Error matrix created externally:
void DifferentialDataStructure::setErrorMatrix( MatrixDataStructure* errm ) {
  if( errm->getBinedges() != binedges ) {
    throw std::runtime_error( "matrix dimension does not match" );
  }
  errorMatrix= errm;
}

// Normalise only bins, not under- or underflow, because binwidth is not defined:
void DifferentialDataStructure::normalise() {
  checkNormalised();
  checkNtotalGTZero();
  if( errorMatrix == 0 ) {
    errorMatrix= new MatrixDataStructure( binedges );
    calculateErrorMatrixWeighted();
  }
  else {
    // Some unfolder already created an error matrix, use Jacobean method
    calculateErrorMatrixJacobean();
    //   cout << "DifferentialDataStructure::normalise: Jacobean method not implemented"
    //	 << endl;
  }    
  for( size_t i= 0; i < points.size(); i++ ) {
    Double_t binw= binedges[i+1]-binedges[i];
    values[i+1]/= (Ntotal*binw);
    errors[i+1]/= (Ntotal*binw);
  }
  setNormalisedTrue();
}

void DifferentialDataStructure::Print() const {
  size_t n= binedges.size();
  cout << "Under- and overflow: " << values[0] << " +/- " << errors[0] 
       << ", " << values[n] << " +/- " << errors[n] << endl;
  for( size_t i= 0; i < n-1; i++ ) {
    cout << binedges[i] << " " << binedges[i+1] << ": " << values[i+1] << " +/- " 
	 << errors[i+1] << endl;
  }
  Double_t sum= 0.0;
  for( size_t i= 0; i < n+1; i++ ) sum+= values[i];
  cout << "Sum of bins incl. under- and overflow: " << sum << endl;
  if( errorMatrix ) {
    cout << "Error matrix:" << endl;
    errorMatrix->Print();
  }
}

// Generalised OPAL PR 404 method for diagonal errors of bin-by-bin corrected
// distribution from calculation of the Jacobean after normalisation.
// This also works for weighted distributions, e.g. y*dsigma/dy or EEC
void DifferentialDataStructure::calculateErrorMatrixWeighted() { 
  size_t nbin= values.size()-2;
  Double_t sumw= 0.0;
  for( size_t ibin= 1; ibin <= nbin; ibin++ ) sumw+= values[ibin];
  for( size_t ibin= 1; ibin <= nbin; ibin++ ) {
    Double_t binwi= binedges[ibin]-binedges[ibin-1];
    for( size_t jbin= 1; jbin <= nbin; jbin++ ) {
      Double_t cov= 0.0;
      for( size_t kbin= 1; kbin <= nbin; kbin++ ) {
        Double_t deltaik= 0.0;
        if( ibin == kbin ) deltaik= 1.0;
        Double_t deltajk= 0.0;
        if( jbin == kbin ) deltajk= 1.0;
        cov+= pow( errors[kbin], 2 )*
          ( sumw*deltaik - values[ibin] )*
          ( sumw*deltajk - values[jbin] );
      }
      cov/= pow( Ntotal*sumw, 2 );
      Double_t binwj= binedges[jbin]-binedges[jbin-1];
      cov/= (binwi*binwj);
      errorMatrix->setElement( ibin, jbin, cov );
    }
  }
  return;
}

void DifferentialDataStructure::calculateErrorMatrixJacobean() {
  size_t nbin= values.size()-2;
  Double_t sumw= 0.0;
  for( size_t ibin= 1; ibin <= nbin; ibin++ ) sumw+= values[ibin];
  MatrixDataStructure Jacobean( binedges );
  for( size_t ibin= 1; ibin <= nbin; ibin++ ) {
    for( size_t jbin= 1; jbin <= nbin; jbin++ ) {
      Double_t deltaij= 0.0;
      if( ibin == jbin ) deltaij= 1.0;
      Jacobean.setElement( ibin, jbin, ( sumw*deltaij - values[ibin] )/(Ntotal*sumw) );
    }
  }
  *errorMatrix= multiply( Jacobean, multiply( *errorMatrix, Jacobean ) );
  for( size_t ibin= 1; ibin <= nbin; ibin++ ) {
    Double_t binwi= binedges[ibin]-binedges[ibin-1];
    for( size_t jbin= 1; jbin <= nbin; jbin++ ) {
      Double_t binwj= binedges[jbin]-binedges[jbin-1];
      Double_t element= errorMatrix->getElement( ibin, jbin );
      errorMatrix->setElement( ibin, jbin, element/(binwi*binwj) );
    }
  }
  return;
}

