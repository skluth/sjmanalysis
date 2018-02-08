
#include "SjmConfigParser.hh"

#include <fstream>
#include <iostream>
#include <exception>
#include <sstream>
#include <boost/algorithm/string.hpp>

SjmConfigParser::SjmConfigParser( int argc, const char* argv[] ) {

  try {

    // Cmd line:
    po::options_description generalOptions( "General" );
    generalOptions.add_options()
      ( "help,h", "Help screen" )
      ( "config", po::value<std::string>(), "Config file" );
    po::positional_options_description posOpts; 
    posOpts.add( "config", 1 );

    // handle help option:
    if( vm.count( "help" ) ) {
      std::cout << generalOptions << std::endl;
    }

    // Parse cmd line to get config filename:
    po::store( po::command_line_parser( argc, argv ).
	       options( generalOptions ).positional( posOpts ).run(),
	       vm );
    if( not vm.count( "config" ) ) {
      throw std::runtime_error( "Config file name not given" );
    }

    // Config file from cmd line"
    po::options_description fileOptions( "File" );
    fileOptions.add_options()
      ( "path", po::value<std::string>()->default_value("."), "Data path" )
      ( "energy", po::value<float>(), "Energy point" )
      ( "maxevt", po::value<int>()->default_value(999999), "Max. number of events" )
      ( "Data.lumi", po::value<float>(), "Data Lumi" )
      ( "Data.files", po::value<std::vector<std::string>>(), 
	"Data file names" )
      ( "SignalMC.xsec", po::value<float>(), "Signal MC xsec" )
      ( "SignalMC.name", po::value<std::string>(), "Signal MC name" )
      ( "SignalMC.files", po::value<std::vector<std::string>>(),
	"Signal MC file names" )
      ( "AltSignalMC.xsec", po::value<float>(), "Alt. signal MC xsec" )
      ( "AltSignalMC.name", po::value<std::string>(), "Alt. signal MC name" )
      ( "AltSignalMC.files", po::value<std::vector<std::string>>(),
	"Alt. signal MC file names" )
      ( "Observables.observable", po::value<std::vector<std::string>>(),
	"Observable names" );
    

    cfgfilename= vm["config"].as<std::string>();
    std::ifstream ifs( cfgfilename );
    if( ifs ) po::store( po::parse_config_file( ifs, fileOptions ), vm );
    else throw std::runtime_error( "file not found" );
    po::notify( vm );

  }
  catch( const po::error & ex ) {
    std::cerr << ex.what() << std::endl;
  }

}

std::string SjmConfigParser::getConfigFileName() const {
  return cfgfilename;
}
std::string SjmConfigParser::getPath() const {
  return getString( "path" );
}
float SjmConfigParser::getEnergy() const {
  return getFloat( "energy" );
}
int SjmConfigParser::getMaxEvent() const {
  return vm["maxevt"].as<int>();
}
float SjmConfigParser::getFloat( const std::string & tag ) const {
  if( vm.count( tag ) ) {
    return vm[tag].as<float>();
  }
  else {
    std::ostringstream txt;
    txt << "SjmConfigParser::getFloat: " << tag << " not found";
    throw std::runtime_error( txt.str() );
  }
}
std::string SjmConfigParser::getString( const std::string & tag ) const {
  if( vm.count( tag ) ) {
    return vm[tag].as<std::string>();
  }
  else {
    std::ostringstream txt;
    txt << "SjmConfigParser::getString: " << tag << " not found";
    throw std::runtime_error( txt.str() );
  }
}
float SjmConfigParser::getDataLumi() const {
  return getFloat( "Data.lumi" );
}
float SjmConfigParser::getSignalMCXsec() const {
  return getFloat( "SignalMC.xsec" );
}
float SjmConfigParser::getAltSignalMCXsec() const {
  return getFloat( "AltSignalMC.xsec" );
}
std::string SjmConfigParser::getSignalMCName() const {
  return getString( "SignalMC.name" );
}
std::string SjmConfigParser::getAltSignalMCName() const {
  return getString( "AltSignalMC.name" );
}
std::vector<std::string> 
SjmConfigParser::getTokens( const std::string& tag ) const {
  if( vm.count( tag ) ) {
    return vm[tag].as<std::vector<std::string>>();
  }
  else {
    std::ostringstream txt;
    txt << "SjmConfigParser::getTokens: " << tag << " not found";
    throw std::runtime_error( txt.str() );
  }
}
std::vector<std::string> 
SjmConfigParser::getFilepath( const std::string & tag ) const {
  std::vector<std::string> files= getTokens( tag );
  std::string path= getPath();
  path+= "/";
  for( std::string & file : files ) {
    file.insert( 0, path );
  }
  return files;
}
std::vector<std::string> SjmConfigParser::getDataFiles() const {
  return getFilepath( "Data.files" );
}
std::vector<std::string> SjmConfigParser::getSignalMCFiles() const {
  return getFilepath( "SignalMC.files" );
}
std::vector<std::string> SjmConfigParser::getAltSignalMCFiles() const {
  return getFilepath( "AltSignalMC.files" );
}
std::vector<std::string> SjmConfigParser::getObservables() const {
  return getTokens( "Observables.observable" );
}



void SjmConfigParser::printTokens( const std::string& txt,
				   const std::vector<std::string> tokens ) const {
  std::cout << txt << std::endl;
  for( size_t i= 0; i < tokens.size(); i++ ) {
    std::cout << " " << i << ": " << tokens[i] << std::endl;
  }
  //std::cout << std::endl;
}

void SjmConfigParser::printConfig() const {

  std::cout << "Configuration: " << getConfigFileName() << std::endl;
  std::cout << "Path: " << getPath() << std::endl;
  std::cout << "Energy: " << getEnergy() << std::endl;
  std::cout << "Max. events: " << getMaxEvent() << std::endl;
  printTokens( "Data files:", getDataFiles() );
  std::cout << "Signal MC name: " << getSignalMCName() << std::endl;
  printTokens( "Signal MC files:", getSignalMCFiles() );
  std::cout << "Alt. signal MC name: " << getAltSignalMCName() << std::endl;
  printTokens( "Alt. Signal MC files:", getAltSignalMCFiles() );
  printTokens( "Observables: ", getObservables() );

}
