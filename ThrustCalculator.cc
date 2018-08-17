
#include "ThrustCalculator.hh"
#include <string>
#include <vector>
#include "Rtypes.h"
#include "TLorentzVector.h"
#include <stdexcept>
#include "NtupleReader.hh"
#include "LEPNtupleReader.hh"

#include <iostream>

// Link PXLIB Fortran
extern "C" {
  void pxlth4_( Int_t*, Int_t*, Float_t*, Float_t*, Float_t*, Int_t* );
};


Double_t ThrustCalculator::getValue( NtupleReader* ntr, 
				     const std::string& reco ) const {  
  static const Int_t maxtrk=500;
  Int_t itkdm= 4;
  Float_t ptrak[maxtrk*itkdm];
  const std::vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( reco );
  Int_t ntrak= vtlv.size();
  for( Int_t itrak= 0; itrak < ntrak; itrak++ ) {
    TLorentzVector tlv= vtlv[itrak];
    Int_t ioff= itkdm*itrak;
    ptrak[ioff]= tlv.Px();
    ptrak[ioff+1]= tlv.Py();
    ptrak[ioff+2]= tlv.Pz();
    ptrak[ioff+3]= tlv.E();    
  }
  Float_t thrval[3];
  Float_t thrvec[3*3];
  Int_t ierr;
  pxlth4_( &ntrak, &itkdm, ptrak, thrval, thrvec, &ierr );
  if( ierr != 0 ) {
    throw std::runtime_error( "ThrustCalculator::getValue: pxlth4 error"+std::to_string( ierr ) );
  }
  return 1.0-thrval[2];
}
