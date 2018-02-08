#ifndef SJMCONFIGPARSER_HH
#define SJMCONFIGPARSER_HH

#include <vector>
#include <string>
#include <boost/program_options.hpp>
namespace po= boost::program_options;

class SjmConfigParser {

  po::variables_map vm;
  std::string cfgfilename;

  float getFloat( const std::string & tag ) const;
  std::string getString( const std::string & tag ) const;
  std::vector<std::string> getTokens( const std::string & tag ) const;
  void printTokens( const std::string& txt,
		    const std::vector<std::string> names ) const;
  std::vector<std::string> getFilepath( const std::string & tag ) const;
  std::string getPath() const;

public:

  SjmConfigParser( int argc, const char* argv[] );
  ~SjmConfigParser() {}

  // po::variables_map getVm() { return vm; }

  std::string getConfigFileName() const;
  float getEnergy() const;
  int getMaxEvent() const;
  float getDataLumi() const;
  float getSignalMCXsec() const;
  float getAltSignalMCXsec() const;
  std::string getSignalMCName() const;
  std::string getAltSignalMCName() const;
  std::vector<std::string> getDataFiles() const;
  std::vector<std::string> getSignalMCFiles() const;
  std::vector<std::string> getAltSignalMCFiles() const;
  void printConfig() const;
  std::vector<std::string> getObservables() const;

};

#endif
