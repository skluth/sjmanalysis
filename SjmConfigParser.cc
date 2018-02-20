
#include "SjmConfigParser.hh"

#include <fstream>
#include <iostream>
#include <exception>
#include <sstream>
#include <boost/algorithm/string.hpp>

SjmConfigParser::SjmConfigParser( int argc, const char* argv[] ) {

  try {

    // Handle cmdline:
    po::options_description cmdlineOptions( "Cmdline" );
    cmdlineOptions.add_options()
      ( "help,h", "Help screen" )
      ( "config", po::value<std::string>(), "Config file" );
    po::positional_options_description posOpts; 
    posOpts.add( "config", 1 );
    po::store( po::command_line_parser( argc, argv ).
	       options( cmdlineOptions ).positional( posOpts ).run(),
	       vm );
    if( vm.count( "help" ) ) {
      std::cout << cmdlineOptions << std::endl;
    }
    if( not vm.count( "config" ) ) {
      throw std::runtime_error( "Config file name not given" );
    }

    // Handle config file from cmd line:
    po::options_description specialOptions( "Special options" );
    specialOptions.add_options()
      ( "General.energy", po::value<std::string>(), "Energy point" )
      ( "General.maxevt", po::value<int>()->default_value(999999), "Max. number of events" )
      ( "General.outfile", po::value<std::string>(), "Output file name" )
      ( "General.sjmGeneralOptions", po::value<std::string>()->
	default_value( "sjmGeneralOptions.cfg" ), "Observable configuration file" )
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
      ( "Analyses.data", po::value<std::vector<std::string>>(), "Data analyses" )
      ( "Analyses.signal", po::value<std::vector<std::string>>(), 
	"Signal MC analyses" )
      ( "Analyses.altsignal", po::value<std::vector<std::string>>(), 
	"Alt. signal MC analyses" );


    cfgfilename= getItem<std::string>( "config" );
    std::ifstream ifs( cfgfilename );
    if( ifs ) po::store( po::parse_config_file( ifs, specialOptions ), vm );
    else throw std::runtime_error( "Config file not found" );
    po::notify( vm );

    // Configuratiom of general options:
    po::options_description generalOptions( "General options" );
    generalOptions.add_options()
      ( "General.url", po::value<std::string>()->default_value( "" ), "Data URL" )
      ( "General.normalise", po::value<bool>()->default_value( true ), "Normalise" )
      ( "General.path", po::value<std::string>()->default_value( "." ), "Data path" )
      ( "Observables.observable", po::value<std::vector<std::string>>(),
	"Observable names" )
      ( "Points.thrust", po::value<std::string>()->default_value( "none" ), 
	"Thrust bin edges" )
      ( "Points.y34cut", po::value<float>()->default_value( 0.0045 ), "y34 cut" )
      ( "Points.y34y23cut", po::value<float>()->default_value( 0.5 ), "y34/y23 cut" )
      ( "Points.a14", po::value<std::string>()->default_value( "none" ), "A14 bin edges" )
      ( "Points.c202", po::value<std::string>()->default_value( "none" ), "C202 bin edges" )
      ( "Points.as", po::value<std::string>()->default_value( "none" ), "AS bin edges" )
      ( "Points.mr", po::value<std::string>()->default_value( "none" ), "MR bin edges" )
      ( "Points.yNMPoints", po::value<std::string>()->default_value( "none" ), "Ynm points" )
      ( "Points.PxEminPoints", po::value<std::string>()->default_value( "none" ), 
	"PXCONE Emin points" )
      ( "Points.PxRPoints", po::value<std::string>()->default_value( "none" ), 
	"PXCONE R points" )
      ( "Points.PxConeEmin", po::value<float>()->default_value( 0.07 ), "PXCONE Emin value" )
      ( "Points.PxConeR", po::value<float>()->default_value( 0.7 ), "PXCONE R value" )
      ( "Points.Donkersycutd", po::value<std::string>()->default_value( "none" ), 
	"Donkers Ycut D values" )
      ( "Points.Donkersycutj", po::value<std::string>()->default_value( "none" ), 
	"Donkers Ycut J values" )
      ( "Points.eminFractionPoints", po::value<std::string>()->default_value( "none" ), 
	"Emin fraction points" )
      ( "Points.RPoints", po::value<std::string>()->default_value( "none" ), "R points" )
      ( "Points.eminFractionValue", po::value<float>()->default_value( 0.06 ), 
	"Emin fraction value" )
      ( "Points.RValue", po::value<float>()->default_value( 0.7 ), "R value" );
    std::ifstream generalOptionsFile( getItem<std::string>( "General.sjmGeneralOptions" ) );
    if( generalOptionsFile ) {
      po::store( po::parse_config_file( generalOptionsFile,
					generalOptions ), vm );
    }
    else {
      throw std::runtime_error( "General options file not found" );
    }
    po::notify( vm );
    
  }
  catch( const po::error & ex ) {
    std::cerr << ex.what() << std::endl;
  }

}

template<typename T>
T SjmConfigParser::getItem( const std::string& tag ) const {
  if( vm.count( tag ) ) {
    return vm[tag].as<T>();
  }
  else {
    std::ostringstream txt;
    txt << "SjmConfigParser::getItem: " << tag << " not found";
    throw std::runtime_error( txt.str() );
  }
}
template bool SjmConfigParser::getItem( const std::string& tag ) const;
template int SjmConfigParser::getItem( const std::string& tag ) const;
template float SjmConfigParser::getItem( const std::string& tag ) const;
template std::string SjmConfigParser::getItem( const std::string& tag ) const;
template std::vector<std::string> SjmConfigParser::getItem( const std::string& tag ) const;
template <>
std::vector<double> SjmConfigParser::getItem( const std::string& tag ) const {
  std::string field= getItem<std::string>( tag );
  if( field == "none" ) {
    std::string txt= "SjmConfigParser::getItem<vector<double>>: field empty for tag " + tag;
    throw std::runtime_error( txt.c_str() );
  }
  std::vector<std::string> tokens;
  boost::split( tokens, field, boost::is_any_of( " " ), boost::token_compress_on );
  std::vector<double> result;
  for( const std::string& token : tokens ) {
    result.push_back( std::stod( token ) );
  }
  return result;
}

std::vector<std::string> 
SjmConfigParser::getFilepath( const std::string & tag ) const {
  std::vector<std::string> files= getItem<std::vector<std::string>>( tag );
  std::string path= getItem<std::string>( "General.path" );
  path+= "/";
  for( std::string & file : files ) {
    file.insert( 0, path );
  }
  return files;
}

std::vector<double>
SjmConfigParser::getPoints( const std::string& tag ) const {
  std::string ptag= "Points."+tag;
  return getItem<std::vector<double>>( ptag );
}

void SjmConfigParser::printConfig() const {
  std::cout << "Configuration: " << getItem<std::string>( "config" ) << std::endl;
  std::cout << "Path: " << getItem<std::string>( "General.path" ) << std::endl;
  std::cout << "Energy: " << getItem<std::string>( "General.energy" ) << std::endl;
  std::cout << "Normalise: ";
  if( getItem<bool>( "General.normalise" ) ) {
    std::cout << "yes";
  }
  else {
    std::cout << "no";
  }
  std::cout << std::endl;    
  std::cout << "Output file: " << getItem<std::string>( "General.outfile" ) << std::endl;
  std::cout << "Max. events: " << getItem<int>( "General.maxevt" ) << std::endl;
  printTokens( "Data files:", getFilepath( "Data.files" ) );
  std::cout << "Signal MC name: " << getItem<std::string>( "SignalMC.name" ) 
	    << std::endl;
  printTokens( "Signal MC files:", getFilepath( "SignalMC.files" ) );
  std::cout << "Alt. signal MC name: " << getItem<std::string>( "AltSignalMC.name" ) 
	    << std::endl;
  printTokens( "Alt. Signal MC files:", getFilepath( "AltSignalMC.files" ) );
  printTokens( "Observables: ", 
	       getItem<std::vector<std::string>>( "Observables.observable" ) );
}

void SjmConfigParser::printTokens( const std::string& txt,
				   const std::vector<std::string> tokens ) const {
  std::cout << txt << std::endl;
  for( size_t i= 0; i < tokens.size(); i++ ) {
    std::cout << " " << i << ": " << tokens[i] << std::endl;
  }
}
