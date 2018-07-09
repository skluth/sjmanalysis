#ifndef TH1DANALYSISOBJECT_HH
#define TH1DANALYSISOBJECT_HH

#include "AnalysisObject.hh"

#include "TString.h"
#include "TVectorD.h"

class TH1D;
class TH2D;

class TH1DAnalysisObject: public AnalysisObject {
public:
  TH1DAnalysisObject( TH1D* h, TH2D* h2d=0 );
  virtual ~TH1DAnalysisObject() {}
  virtual TVectorD getErrors( const TString & opt="" );
  virtual TString getPointStr( Int_t i );
  virtual TString getPointLabel();
  virtual TVectorD getPointsCenter();
private:
  TH1D* hist;
  TH2D* hist2d;
};

#endif
