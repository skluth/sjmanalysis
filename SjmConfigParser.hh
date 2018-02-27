#ifndef SJMCONFIGPARSER_HH
#define SJMCONFIGPARSER_HH

#include <vector>
#include <string>
#include <boost/program_options.hpp>
namespace po= boost::program_options;

class SjmConfigParser {

  po::options_description cmdlineOptions;
  po::options_description specialOptions;
  po::options_description generalOptions;

  po::variables_map vm;
  std::string cfgfilename;

  void setDefault( const std::string& key,
		   const std::string& value );
  void printOption(  boost::shared_ptr<po::option_description> option ) const;
  
public:

  SjmConfigParser( int argc, const char* argv[] );
  ~SjmConfigParser() {}

  template <typename T> T getItem( const std::string& tag ) const;
  std::vector<std::string> getFilepath( const std::string & tag ) const;
  std::vector<double> getPoints( const std::string& tag ) const;
  void printConfig() const;

};

#endif
