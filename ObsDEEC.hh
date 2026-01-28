#ifndef OBSDEEC_HH
#define OBSDEEC_HH

#include "Rtypes.h"
#include "ObsEEC.hh"
#include "Analysis.hh"
#include <vector>

class TLorentzVector;

class ObsDEEC : public ObsEEC {

public:

  ObsDEEC( const std::string & name,
	   const std::vector<Double_t> & bins,
	   const std::vector<Analysis> & variations,
	   const bool scOpt=true,
	   const bool lprint=true );
  virtual ~ObsDEEC();
  
private:

  virtual void calcWeight( const TLorentzVector&, const TLorentzVector&,
			   Double_t&, Double_t& );
  
};

#endif
