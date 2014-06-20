#ifndef JETRATEDATASTRUCTURE_HH
#define JETRATEDATASTRUCTURE_HH

#include "Rtypes.h"
#include "DataStructure.hh"

#include <vector>
using std::vector;

class JetrateDataStructure : public DataStructure {

public:

  JetrateDataStructure( const vector<Double_t>&, Int_t );
  JetrateDataStructure() {}
  ~JetrateDataStructure() {}

  void fill( const vector<Double_t>& NJets );
  void normalise();

  const vector<Double_t>& getPoints() const { return points; }
  void print();
  DataStructure* clone();

private:

  Int_t Jetrate;

};

#endif
