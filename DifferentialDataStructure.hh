#ifndef DIFFERENTIALDATASTRUCTURE_HH
#define DIFFERENTIALDATASTRUCTURE_HH

#include "Rtypes.h"
#include "DataStructure.hh"

#include <vector>
#include <string>


class DifferentialDataStructure : public DataStructure {

public:

  DifferentialDataStructure( const std::vector<Double_t> & p );
  DifferentialDataStructure() {}
  virtual ~DifferentialDataStructure() {}

  void fill( Double_t value, Double_t weight=1.0 );
  void normalise();

  std::vector<Double_t> getBinedges() const { return binedges; }
  void Print() const;
  void printBinedges();
  DataStructure* clone() const;
  virtual void setErrorMatrix();


private:

  void calculateErrorMatrixWeighted();

  std::vector<Double_t> binedges;

};



#endif
