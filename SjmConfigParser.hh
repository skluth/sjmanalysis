#ifndef SJMCONFIGPARSER_HH
#define SJMCONFIGPARSER_HH

#include <vector>
#include <string>
#include <boost/program_options.hpp>
namespace po= boost::program_options;

class SjmConfigParser {

  po::variables_map vm;
  std::string cfgfilename;

  void printTokens( const std::string& txt,
		    const std::vector<std::string> names ) const;
  //std::string getPath() const;

public:

  SjmConfigParser( int argc, const char* argv[] );
  ~SjmConfigParser() {}

  template <typename T> T getItem( const std::string& tag ) const;
  std::vector<std::string> getFilepath( const std::string & tag ) const;
  std::vector<double> getPoints( const std::string& tag ) const;
  void printConfig() const;

};

#endif
