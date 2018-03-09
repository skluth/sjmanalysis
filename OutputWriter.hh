#ifndef OUTPUTWRITER_HH
#define OUTPUTWRITER_HH

#include <vector>
#include <string>

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

private:

  void writeJetrate( const JetrateDataStructure*,
		     const std::string & );
  void writeDifferentialDistribution( const DifferentialDataStructure*,
				      const std::string & );
  void writeMatrix( MatrixDataStructure*, const std::string & );

  TFile* outputfile;

};



#endif
