#ifndef ANALYSISOBJECT_HH
#define ANALYSISOBJECT_HH

#include "TVectorD.h"
#include "TString.h"

class AnalysisObject {
public:
  AnalysisObject( Int_t n ) : points(n), values(n), errors(n) {}
  virtual ~AnalysisObject() {}
  virtual TVectorD getPoints() { return points; }
  virtual TVectorD getValues() { return values; }
  virtual TVectorD getErrors() { return errors; }
  virtual TString getPointStr( Int_t ) = 0;
  virtual TVectorD getPointsCenter() = 0;
protected:
  TVectorD points;
  TVectorD values;
  TVectorD errors;
};

#endif
