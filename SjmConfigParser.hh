#ifndef SJMCONFIGPARSER_HH
#define SJMCONFIGPARSER_HH

#include <vector>
#include <string>
#include <map>

class SjmConfigParser {

  std::map< std::string, std::string > descriptionMap;
  std::map< std::string, std::vector<std::string> > valuesMap;

  void declareOptions( const std::string & filename,
		       const std::map< std::string, std::string > & descriptions );

public:

  SjmConfigParser( int argc, const char* argv[] );
  ~SjmConfigParser() {}

  template <typename T> T getItem( const std::string& tag ) const;
  std::vector<std::string> getFilepath( const std::string & tag ) const;
  std::vector<double> getPoints( const std::string& tag ) const;
  void printConfig() const;
 
};

#endif
