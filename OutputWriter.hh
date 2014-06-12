#ifndef OUTPUTWRITER_HH
#define OUTPUTWRITER_HH

#include <vector>
using std::vector;
#include <string>
using std::string;
class TFile;

class Observable;

class OutputWriter {

public:

  OutputWriter( const string& filename );
  ~OutputWriter();

  void write( const vector<Observable*>& );

private:

  TFile* outputfile;

};



#endif
