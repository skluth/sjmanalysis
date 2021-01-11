
#include "SjmConfigParser.hh"
#include "AnalysisProcessor.hh"
#include <iostream>

int main( int argc, const char *argv[] ) {

  SjmConfigParser sjmcp( argc, argv );
  sjmcp.printConfig();
  AnalysisProcessor ap( sjmcp );
  ap.LEPAnalysis();
  return 0;

}

