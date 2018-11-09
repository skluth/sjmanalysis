#ifndef TGEANALYSISOBJECT_HH
#define TGEANALYSISOBJECT_HH

#include "AnalysisObject.hh"
#include "TString.h"
#include "TVectorD.h"

class TGraphErrors;

class TGEAnalysisObject: public AnalysisObject {
public:
  TGEAnalysisObject( TGraphErrors* t );
  virtual ~TGEAnalysisObject() {}
  virtual TString getPointStr( Int_t i, Int_t width=5, Int_t prec=2 );
  virtual TString getPointLabel( Int_t width=5 );
  virtual TVectorD getPointsCenter();
};

#endif
