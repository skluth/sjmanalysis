#ifndef PXTHRUSTCALCULATOR_HH
#define PXTHRUSTCALCULATOR_HH

#include "DifferentialCalculator.hh"

class NtupleReader;
#include <string>

class PxThrustCalculator: public DifferentialCalculator {
  static const Int_t maxtrk= 500;
  static const Int_t itkdm= 4;
  // need these non-const in getValue for Fortran interface
  mutable Int_t itkdmpx= itkdm;
  mutable Float_t ptrak[maxtrk*itkdm];
public:
  PxThrustCalculator() {}
  ~PxThrustCalculator() {}
  Double_t getValue( NtupleReader*, const std::string & ) const;
};

#endif
