#ifndef JETRATEDATASTRUCTURE_HH
#define JETRATEDATASTRUCTURE_HH

#include "Rtypes.h"
#include "DataStructure.hh"

#include <vector>
using std::vector;

class JetrateDataStructure : public DataStructure {

public:

  JetrateDataStructure( const vector<Double_t>& p );
  JetrateDataStructure() {}
  ~JetrateDataStructure() {}

  void fill( const vector<Double_t>& NJets, Int_t Jetrate );
  void normalise();

  const vector<Double_t>& getPoints() const { return points; }
  vector<Double_t> getValues() const { return values; }
  vector<Double_t> getErrors() const { return errors; }
  void setValues( const vector<Double_t>& valuesin ) { values= valuesin; }
  void setErrors( const vector<Double_t>& errorsin ) { errors= errorsin; }
  void print();
  DataStructure* clone();

private:

  vector<Double_t> points;
  vector<Double_t> values;
  vector<Double_t> errors;

};

#endif
