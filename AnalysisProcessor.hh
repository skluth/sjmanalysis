#ifndef ANALYSISPROCESSOR_HH
#define ANALYSISPROCESSOR_HH

#include "Analysis.hh"
#include "SjmConfigParser.hh"

#include <string>
#include <vector>
#include <map>
#include "Rtypes.h"

class Observable;
class FilledObservable;
class LEPNtupleReader;
class NtupleReader;

class AnalysisProcessor {

  SjmConfigParser sjmConfigs;
  int maxevt;

  Int_t processAnalyses( const std::vector<Analysis> & analyses,
  			 const std::vector<Observable*> & vobs,
  			 const std::string & filename );
  Int_t processAnalysesNtr( const std::vector<Analysis> & analyses,
			    const std::vector<Observable*> & vobs,
			    NtupleReader* ntr );

  void processUnfolding( const std::vector<Analysis>& measuredAnalyses, 
			 const std::string& unfoldsource,
			 const std::vector<FilledObservable*>& vfobs );

  std::vector<FilledObservable*> getFilled( const std::vector<Observable*>& vobs );

  std::vector<Analysis> fillAnalyses( const std::string& tag );

  NtupleReader* createNtupleReader( const std::string& filename );

  std::vector<Analysis>
  subtractBackground( const std::vector<FilledObservable*> & vfobs,
		      const std::vector<Analysis> & measuredAnalyses, 
		      const std::map<std::string,Double_t> & eventCounts );

public:

  AnalysisProcessor( const SjmConfigParser& sjmcp );

  ~AnalysisProcessor() {}

  void LEPAnalysis();
  void MCAnalysis();

};

#endif
