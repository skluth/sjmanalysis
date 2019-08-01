#ifndef OUTPUTWRITER_HH
#define OUTPUTWRITER_HH

#include <vector>
#include <string>
#include <map>

class TFile;
class FilledObservable;
class JetrateDataStructure;
class DifferentialDataStructure;
class MatrixDataStructure;

class OutputWriter {

public:

  OutputWriter( const std::string & filename );
  ~OutputWriter();

  void write( const std::vector<FilledObservable*> & );
  
  void writeCutflowCounters( const std::map<std::string,std::map<std::string,int>> & );

private:

  void writeJetrate( const JetrateDataStructure*,
		     const std::string & );
  void writeDifferentialDistribution( const DifferentialDataStructure*,
				      const std::string & );
  void writeMatrix( MatrixDataStructure*, const std::string &, const std::string & );

  TFile* outputfile;

};



#endif
