#ifndef DATASTRUCTURE_HH
#define DATASTRUCTURE_HH

#include "Rtypes.h"

#include <vector>
using std::vector;

class DataStructure {

public:

  DataStructure() : Ntotal(0) {}
  virtual ~DataStructure() {}

  virtual const vector<Double_t>& getPoints() const { return points; }
  virtual const vector<Double_t>& getValues() const { return values; }
  virtual const vector<Double_t>& getErrors() const { return errors; }
  virtual void setValues( const vector<Double_t>& v ) { values=v; }
  virtual void setErrors( const vector<Double_t>& e ) { errors=e; }
  void setNEvents( Double_t nevents ) { Ntotal= nevents; }
  Double_t getNEvents() const { return Ntotal; }
  virtual void normalise() = 0;
  virtual void print() = 0;
  virtual DataStructure* clone() = 0;

protected:

  Double_t Ntotal;
  vector<Double_t> points;
  vector<Double_t> values;
  vector<Double_t> errors;

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
