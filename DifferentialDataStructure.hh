#ifndef DIFFERENTIALDATASTRUCTURE_HH
#define DIFFERENTIALDATASTRUCTURE_HH

#include "Rtypes.h"
#include "DataStructure.hh"

#include <vector>
using std::vector;

class DifferentialDataStructure : public DataStructure {

public:

  DifferentialDataStructure( const vector<Double_t>& p );
  DifferentialDataStructure() {}
  ~DifferentialDataStructure() {}

  void fill( Double_t value, Double_t weight=1.0 );
  void normalise();

  const vector<Double_t>& getPoints() const { return points; }
  vector<Double_t> getValues() const;
  vector<Double_t> getErrors() const;
  const vector<Double_t>& getBinedges() const { return binedges; }
  void setValues( const vector<Double_t>& );
  void setErrors( const vector<Double_t>& );
  void print();
  DataStructure* clone();

private:

  vector<Double_t> binedges;
  vector<Double_t> points;
  vector<Double_t> values;
  vector<Double_t> errors;

};



#endif
