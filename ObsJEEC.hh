#ifndef OBSJEEC_HH
#define OBSJEEC_HH

#include "Rtypes.h"
#include "ObsEEC.hh"
#include "Analysis.hh"
#include <vector>

class TLorentzVector;

class ObsJEEC : public ObsEEC {

public:

  ObsJEEC( const std::string & name,
	   const std::vector<Double_t> & bins,
	   const std::vector<Analysis> & variations,
	   const bool scOpt=true,
	   const bool lprint=true );
  virtual ~ObsJEEC();
  
private:

  virtual void calcWeight( const TLorentzVector&, const TLorentzVector&,
			   Double_t&, Double_t& );
  
};

#endif
