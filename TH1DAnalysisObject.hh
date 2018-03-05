#ifndef TH1DANALYSISOBJECT_HH
#define TH1DANALYSISOBJECT_HH

#include "AnalysisObject.hh"

#include "TString.h"
#include "TVectorD.h"

class TH1D;

class TH1DAnalysisObject: public AnalysisObject {
public:
  TH1DAnalysisObject( TH1D* h );
  virtual ~TH1DAnalysisObject() {}
  virtual TString getPointStr( Int_t i );
  virtual TVectorD getPointsCenter();
private:
  TH1D* hist;
};

#endif
