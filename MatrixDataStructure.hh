#ifndef MATRIXDATASTRUCTURE_HH
#define MATRIXDATASTRUCTURE_HH

#include "Rtypes.h"
#include <vector>

class MatrixDataStructure {

public:

  MatrixDataStructure( const std::vector<Double_t> & bins );
  MatrixDataStructure( const MatrixDataStructure & m );
  virtual ~MatrixDataStructure();
  MatrixDataStructure & operator= ( const MatrixDataStructure & other );

  MatrixDataStructure* clone() const;
  void fill( Double_t, Double_t );
  Double_t getElement( size_t irow, size_t icol ) const;
  void setElement( size_t irow, size_t icol, Double_t value );
  void Print() const;
  Double_t getNEvents() const { return Ntotal; }
  std::vector<Double_t> getBinedges() const { return binedges; } 
  void normaliseColumns();
  size_t getNdim() const { return ndim; }

private:

  void checkIndices( size_t irow, size_t icol ) const;

  std::vector<Double_t> binedges;
  Int_t Ntotal;
  size_t ndim;
  Double_t* array;

};

std::vector<Double_t> multiply( const MatrixDataStructure & m,
				const std::vector<Double_t> & v );
MatrixDataStructure multiply( const MatrixDataStructure & mlhs,
			      const MatrixDataStructure & mrhs );

MatrixDataStructure* applyEfficiency( const MatrixDataStructure & m,
				      const std::vector<Double_t> & v );

MatrixDataStructure* similarityVector( const MatrixDataStructure & m,
				       const std::vector<Double_t> & v );



#endif
