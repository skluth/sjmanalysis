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
  virtual TString getPointStr( Int_t i );
  virtual TString getPointLabel();
  virtual TVectorD getPointsCenter();
private:
  TGraphErrors* tge;
};

#endif
