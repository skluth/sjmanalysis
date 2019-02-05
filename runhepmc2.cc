
#include "AnalysisProcessor.hh"

int main( int argc, const char *argv[] ) {
  SjmConfigParser sjmcp( argc, argv );
  sjmcp.printConfig();
  AnalysisProcessor ap( sjmcp );
  ap.MCAnalysis();
  return 0;
}
