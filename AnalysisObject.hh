#ifndef ANALYSISOBJECT_HH
#define ANALYSISOBJECT_HH

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TString.h"
#include <iostream>

class AnalysisObject {
public:
  AnalysisObject( UInt_t n ) : points(n), values(n), errors(n) {}
  virtual ~AnalysisObject() {}
  virtual TVectorD getPoints() { return points; }
  virtual TVectorD getValues() { return values; }
  virtual TVectorD getErrors( const TString & opt="" ) {
    TVectorD result= errors;
    if( opt.Index( "m" ) >= 0 ) {
      if( errorMatrix.GetNoElements() > 0 ) {
	for( Int_t i= 0; i < errorMatrix.GetNrows(); i++ ) {
	  result[i]= sqrt( errorMatrix[i][i] );
	}
      }
      else {
	std::cout << "AnalysisObject::getErrors: error matrix empty" << std::endl;
      }
    }
    return result;
  }
  virtual TMatrixD getErrorMatrix() { return errorMatrix; }
  virtual TString getPointStr( Int_t, Int_t, Int_t ) = 0;
  virtual TString getPointLabel( Int_t ) = 0;
  virtual TVectorD getPointsCenter() = 0;
protected:
  TVectorD points;
  TVectorD values;
  TVectorD errors;
  TMatrixD errorMatrix;
};

#endif
