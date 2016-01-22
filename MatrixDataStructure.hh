#ifndef MATRIXDATASTRUCTURE_HH
#define MATRIXDATASTRUCTURE_HH

#include "Rtypes.h"
#include <vector>
using std::vector;

class MatrixDataStructure {

public:

  MatrixDataStructure( const vector<Double_t>& bins );
  ~MatrixDataStructure() {}

  void fill( Double_t, Double_t );
  Double_t getElement( size_t irow, size_t icol ) const;
  void print() const;
  Double_t getNEvents() const { return Ntotal; }
  const vector<Double_t>& getBinedges() const { return binedges; } 

private:

  vector<Double_t> binedges;
  Int_t Ntotal;
  size_t ndim;
  Double_t* array;

};

#endif
