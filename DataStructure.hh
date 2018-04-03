#ifndef DATASTRUCTURE_HH
#define DATASTRUCTURE_HH

#include "Rtypes.h"

#include <vector>

class MatrixDataStructure;

class DataStructure {

public:

  DataStructure() : Ntotal(0.0), errorMatrix(0), normalised(false) {}
  virtual ~DataStructure() {}

  virtual const std::vector<Double_t>& getPoints() const { return points; }
  virtual const std::vector<Double_t>& getValues() const { return values; }
  virtual const std::vector<Double_t>& getErrors() const { return errors; }
  virtual void setValues( const std::vector<Double_t>& v ) { values=v; }
  virtual void setErrors( const std::vector<Double_t>& e ) { errors=e; }
  void setNEvents( Double_t nevents ) { Ntotal= nevents; }
  Double_t getNEvents() const { return Ntotal; }
  void checkNormalised();
  void setNormalisedTrue() { normalised= true; }
  bool getNormalised() const { return normalised; }
  void checkNtotalGTZero();
  virtual void normalise() = 0;
  virtual void Print() const = 0;
  virtual DataStructure* clone() const = 0;
  virtual void setErrorMatrix( MatrixDataStructure* ) = 0;
  MatrixDataStructure* getErrorMatrix() const;
  
protected:

  Double_t Ntotal;
  std::vector<Double_t> points;
  std::vector<Double_t> values;
  std::vector<Double_t> errors;
  MatrixDataStructure* errorMatrix;

  
private:

  bool normalised;

};

#endif
