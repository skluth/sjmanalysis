#ifndef DATASTRUCTURE_HH
#define DATASTRUCTURE_HH

#include "Rtypes.h"

#include <vector>
using std::vector;

class DataStructure {

public:

  DataStructure() : Ntotal(0) {}
  ~DataStructure() {}

  virtual void normalise() = 0;

  virtual const vector<Double_t>& getPoints() const = 0;
  virtual vector<Double_t> getValues() const = 0;
  virtual vector<Double_t> getErrors() const = 0;
  virtual void setValues( const vector<Double_t>& ) = 0;
  virtual void setErrors( const vector<Double_t>& ) = 0;
  void setNEvents( Double_t nevents ) { Ntotal= nevents; }
  Double_t getNEvents() const { return Ntotal; }
  virtual void print() = 0;
  virtual DataStructure* clone() = 0;

protected:

  Double_t Ntotal;

};

template <typename T>
vector<T> divideVectors( const vector<T>& lhs, const vector<T>& rhs ) { 
  size_t n= lhs.size();
  vector<T> result( n );
  for( size_t i= 0; i < n; i++ ) {
    if( rhs[i] != 0.0 ) result[i]= lhs[i]/rhs[i];
    else result[i]= 0.0;
  }
  return result;
}
template <typename T>
vector<T> multiplyVectors( const vector<T>& lhs, const vector<T>& rhs ) {
  size_t n= lhs.size();
  vector<T> result( n );
  for( size_t i= 0; i < n; i++ ) result[i]= lhs[i]*rhs[i];
  return result;
}
template <typename T>
vector<T> subtractVectors( const vector<T>& lhs, const vector<T>& rhs ) {
  size_t n= lhs.size();
  vector<T> result( n );
  for( size_t i= 0; i < n; i++ ) result[i]= lhs[i]-rhs[i];
  return result;
}

#endif
