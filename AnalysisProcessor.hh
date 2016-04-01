#ifndef ANALYSISPROCESSOR_HH
#define ANALYSISPROCESSOR_HH

#include "Analysis.hh"
#include <string>
#include <vector>

class Observable;
class FilledObservable;
class SjmConfigParser;

class AnalysisProcessor {

  std::vector<std::string> datafilenames;
  std::vector<std::string> pyfilenames;
  std::vector<std::string> hwfilenames;
  int maxevt;
  std::vector<std::string> obsnames;

  void processAnalyses( const std::vector<Analysis>& analyses,
			const std::vector<Observable*>& vobs,
			const std::string& filename );

  void processUnfolding( const std::vector<Analysis>& measuredAnalyses, 
			 const std::string& unfoldsource,
			 const std::vector<FilledObservable*>& vobs );

  std::vector<FilledObservable*> getFilled( const std::vector<Observable*>& vobs );

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
