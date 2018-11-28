
#include "PxThrustCalculator.hh"
#include <string>
#include <vector>
#include "Rtypes.h"
#include "TLorentzVector.h"
#include <stdexcept>
#include "NtupleReader.hh"

// Link PXLIB thrust calculation
extern "C" {
  void pxlth4_( Int_t*, Int_t*, Float_t*, Float_t*, Float_t*, Int_t* );
};

Double_t PxThrustCalculator::getValue( NtupleReader* ntr, 
				       const std::string& reco ) const {
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
  pxlth4_( &ntrak, &itkdmpx, ptrak, thrval, thrvec, &ierr );
  if( ierr != 0 ) {
    throw std::runtime_error( "PxThrustCalculator::getValue: pxlth4 error "+std::to_string( ierr ) );
  }
  return 1.0-thrval[2];
}

