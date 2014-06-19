#ifndef OUTPUTWRITER_HH
#define OUTPUTWRITER_HH

#include <vector>
using std::vector;
#include <string>
using std::string;
class TFile;

class FilledObservable;
class JetrateDataStructure;
class DifferentialDataStructure;
class MatrixDataStructure;

class OutputWriter {

public:

  OutputWriter( const string& filename );
  ~OutputWriter();

  void write( const vector<FilledObservable*>& );

private:

  void writeJetrate( JetrateDataStructure*, const string& );
  void writeDifferentialDistribution( DifferentialDataStructure*, const string& );
  void writeMatrix( MatrixDataStructure*, const string& );

  TFile* outputfile;

};



#endif
