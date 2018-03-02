#ifndef ANALYSISPROCESSOR_HH
#define ANALYSISPROCESSOR_HH

#include "Analysis.hh"
#include "SjmConfigParser.hh"

#include <string>
#include <vector>
#include "Rtypes.h"

class Observable;
class FilledObservable;
class NtupleReader;

class AnalysisProcessor {

  SjmConfigParser sjmConfigs;
  int maxevt;

  Int_t processAnalyses( const std::vector<Analysis>& analyses,
			 const std::vector<Observable*>& vobs,
			 const std::string& filename );

  void processUnfolding( const std::vector<Analysis>& measuredAnalyses, 
			 const std::string& unfoldsource,
			 const std::vector<FilledObservable*>& vobs );

  std::vector<FilledObservable*> getFilled( const std::vector<Observable*>& vobs );

  std::vector<Analysis> fillAnalyses( const std::string& tag );

  NtupleReader* createNtupleReader( const string& filename );

public:

  // AnalysisProcessor() : datafilename( "da91_96_200.root" ),
  // 			pyfilename( "mc5025_1_200.root" ), 
  //                       hwfilename( "mc12406_1_200.root" ),
  // 		        maxevt( 1000 ) {}
  AnalysisProcessor( const SjmConfigParser& sjmcp );
  ~AnalysisProcessor() {}

  void LEP1Analysis();

};

#endif
