
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

// Do the calculation, store results in thrval and thrvec
void PxThrustCalculator::calculate( const std::vector<TLorentzVector> & vtlv ) const {
  Int_t ntrak= vtlv.size();
  for( Int_t itrak= 0; itrak < ntrak; itrak++ ) {
    TLorentzVector tlv= vtlv[itrak];
    Int_t ioff= itkdm*itrak;
    ptrak[ioff]= tlv.Px();
    ptrak[ioff+1]= tlv.Py();
    ptrak[ioff+2]= tlv.Pz();
    ptrak[ioff+3]= tlv.E();    
  }
  Int_t ierr;
  pxlth4_( &ntrak, &itkdmpx, ptrak, thrval, thrvec, &ierr );
  if( ierr != 0 ) {
    throw std::runtime_error( "PxThrustCalculator::calculate: pxlth4 error " +
			      std::to_string( ierr ) );
  }
  return;
}

// Return results
Double_t PxThrustCalculator::getValue( NtupleReader* ntr, 
				       const std::string & reco ) const {
  const std::vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( reco );
  calculate( vtlv );
  return 1.0-thrval[2];
}
// TVector3 PxThrustCalculator::getTVector( NtupleReader* ntr,
// 					 const std::string & reco ) const {
//   const std::vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( reco );
//   calculate( vtlv );
//   // last three elements of thrvec[3*3] array are Thrust principal axis
//   TVector3 tvec( thrvec[6], thrvec[7], thrvec[8] );
//   return tvec;
// }
TVector3 PxThrustCalculator::getTVector( const std::vector<TLorentzVector> & vtlv ) const {
  calculate( vtlv );
  // last three elements of thrvec[3*3] array are Thrust principal axis
  TVector3 tvec( thrvec[6], thrvec[7], thrvec[8] );
  return tvec;
}
