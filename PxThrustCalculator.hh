#ifndef PXTHRUSTCALCULATOR_HH
#define PXTHRUSTCALCULATOR_HH

#include "DifferentialCalculator.hh"

class NtupleReader;
#include <string>
#include <vector>

#include "TLorentzVector.h"
#include "TVector3.h"

class PxThrustCalculator: public DifferentialCalculator {
  static const Int_t maxtrk= 500;
  static const Int_t itkdm= 4;
  // need these non-const in calculate for Fortran interface
  mutable Int_t itkdmpx= itkdm;
  mutable Float_t ptrak[maxtrk*itkdm];
  mutable Float_t thrval[3];
  mutable Float_t thrvec[3*3];
  void calculate( const std::vector<TLorentzVector> & ) const;
public:
  PxThrustCalculator() {}
  ~PxThrustCalculator() {}
  Double_t getValue( NtupleReader*, const std::string & ) const;
  // TVector3 getTVector( NtupleReader*, const std::string & ) const;
  TVector3 getTVector( const std::vector<TLorentzVector> & ) const;
};

#endif
