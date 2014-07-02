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

  const vector<Double_t>& getBinedges() const { return binedges; }
  void print();
  void printBinedges();
  DataStructure* clone();

private:

  vector<Double_t> binedges;

};



#endif
