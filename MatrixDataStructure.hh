#ifndef MATRIXDATASTRUCTURE_HH
#define MATRIXDATASTRUCTURE_HH

#include "Rtypes.h"
#include <vector>

class MatrixDataStructure {

public:

  MatrixDataStructure( const std::vector<Double_t> & bins );
  ~MatrixDataStructure() {}

  MatrixDataStructure* clone() const;
  void fill( Double_t, Double_t );
  Double_t getElement( size_t irow, size_t icol ) const;
  void setElement( size_t irow, size_t icol, Double_t value );
  void Print() const;
  Double_t getNEvents() const { return Ntotal; }
  std::vector<Double_t> getBinedges() const { return binedges; } 

private:

  void checkIndices( size_t irow, size_t icol ) const;

  std::vector<Double_t> binedges;
  Int_t Ntotal;
  size_t ndim;
  Double_t* array;

};

#endif
