#ifndef JETRATEDATASTRUCTURE_HH
#define JETRATEDATASTRUCTURE_HH

#include "Rtypes.h"
#include "DataStructure.hh"

#include <vector>
#include <string>


class JetrateDataStructure : public DataStructure {

public:

  JetrateDataStructure( const std::vector<Double_t>&, Int_t, const std::string & opt="?" );
  JetrateDataStructure() {}
  virtual ~JetrateDataStructure() {}

  void fill( const std::vector<Double_t>& NJets );
  void normalise();

  const std::vector<Double_t>& getPoints() const { return points; }
  void Print() const;
  DataStructure* clone() const;
  virtual void setErrorMatrix( MatrixDataStructure* );

private:

  Int_t Jetrate;
  std::string inclopt;
  
};

#endif
