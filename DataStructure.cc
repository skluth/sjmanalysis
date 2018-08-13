
#include "DataStructure.hh"

#include <sstream>
using std::ostringstream;
#include <stdexcept>

void DataStructure::checkNormalised() {
  if( normalised ) {
    ostringstream txt;
    txt << "Already normalised";
    throw std::runtime_error( txt.str() );
  }
}

void DataStructure::checkNtotalGTZero() {
  if( Ntotal <= 0.0 ) {
    ostringstream txt;
    txt << "Ntotal <= 0.0: " << Ntotal;
    throw std::runtime_error( txt.str() );
  }
}

MatrixDataStructure* DataStructure::getErrorMatrix() const {
  return errorMatrix;
}
