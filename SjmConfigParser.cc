
#include "SjmConfigParser.hh"

#include <fstream>
#include <iostream>
#include <exception>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>

using std::cout;
using std::endl;

typedef std::vector<std::string> VS;

template <>
std::string SjmConfigParser::getItem( const std::string& tag ) const;

SjmConfigParser::SjmConfigParser( int argc, const char* argv[] ) :
  cmdlineOptions( "Cmdline" ),
  specialOptions( "Special options" ),
  generalOptions( "General options" ) {

  try {

    // Handle cmdline:
    cmdlineOptions.add_options()
      ( "help,h", "Help screen" )
      ( "config", po::value<VS>(), "Config file" );
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
    specialOptions.add_options()

      ( "General.test", po::value<VS>(), "test field" )

      ( "General.energy", po::value<VS>(), "Energy point" )
      ( "General.maxevt", po::value<VS>(), "Max. number of events" )
      ( "General.outfile", po::value<VS>(), "Output file name" )
      ( "General.sjmGeneralOptions", po::value<VS>(), "Observable configuration file" )
      ( "Data.lumi", po::value<VS>(), "Data Lumi" )
      ( "Data.files", po::value<VS>(), "Data file names" )
      ( "Signal.xsec", po::value<VS>(), "Signal MC xsec" )
      ( "Signal.name", po::value<VS>(), "Signal MC name" )
      ( "Signal.files", po::value<VS>(), "Signal MC file names" )
      ( "AltSignal.xsec", po::value<VS>(), "Alt. signal MC xsec" )
      ( "AltSignal.name", po::value<VS>(), "Alt. signal MC name" )
      ( "AltSignal.files", po::value<VS>(), "Alt. signal MC file names" )

      ( "BkgWWllqq.xsec", po::value<VS>(), "Bkg WW->llqq MC xsec" )
      ( "BkgWWllqq.name", po::value<VS>(), "Bkg WW->llqq MC name" )
      ( "BkgWWllqq.files", po::value<VS>(), "Bkg WW->llqq MC file names" )
      
      ( "BkgWWqqqq.xsec", po::value<VS>(), "Bkg WW->qqqq MC xsec" )
      ( "BkgWWqqqq.name", po::value<VS>(), "Bkg WW->qqqq MC name" )
      ( "BkgWWqqqq.files", po::value<VS>(), "Bkg WW->qqqq MC file names" )

      ( "BkgWWeeqq.xsec", po::value<VS>(), "Bkg WW->eeqq MC xsec" )
      ( "BkgWWeeqq.name", po::value<VS>(), "Bkg WW->eeqq MC name" )
      ( "BkgWWeeqq.files", po::value<VS>(), "Bkg WW->eeqq MC file names" )      

      ( "Analyses.data", po::value<VS>(), "Data analyses" )
      ( "Analyses.signal", po::value<VS>(), "Signal MC analyses" )
      ( "Analyses.altsignal", po::value<VS>(), "Alt. signal MC analyses" )

      ( "Analyses.bkgllqq", po::value<VS>(), "Bkg WW->llqq MC analyses" )
      ( "Analyses.bkgqqqq", po::value<VS>(), "Bkg WW->qqqq MC analyses" )
      ( "Analyses.bkgeeqq", po::value<VS>(), "Bkg WW->eeqq MC analyses" );

    cfgfilename= getItem<std::string>( "config" );
    std::ifstream ifs( cfgfilename );
    if( ifs ) po::store( po::parse_config_file( ifs, specialOptions ), vm );
    else throw std::runtime_error( "Config file not found" );
    setDefault( "General.test", "bla" );
    setDefault( "General.maxevt", "999999" );
    setDefault( "General.sjmGeneralOptions", "sjmGeneralOptions.cfg" );
    po::notify( vm );
    
    // Configuratiom of general options:
    generalOptions.add_options()
      ( "General.url", po::value<VS>(), "Data URL" )
      ( "General.normalise", po::value<VS>(), "Normalise" )
      ( "General.path", po::value<VS>(), "Data path" )
      ( "Observables.observable", po::value<VS>(), "Observable names" )

      ( "Points.thrust", po::value<VS>(), "Thrust bin edges" )      
      ( "Points.EEC", po::value<VS>(), "EEC bin edges" )

      ( "Points.y34cut", po::value<VS>(), "y34 cut" )
      ( "Points.y34y23cut", po::value<VS>(), "y34/y23 cut" )
      ( "Points.a14", po::value<VS>(), "A14 bin edges" )
      ( "Points.c202", po::value<VS>(), "C202 bin edges" )
      ( "Points.as", po::value<VS>(), "AS bin edges" )
      ( "Points.mr", po::value<VS>(), "MR bin edges" )

      ( "Points.yNMPoints", po::value<VS>(), "Ynm points" )
      ( "Points.PxEminPoints", po::value<VS>(), "PXCONE Emin points" )
      ( "Points.PxRPoints", po::value<VS>(), "PXCONE R points" )
      ( "Points.PxConeEmin", po::value<VS>(), "PXCONE Emin value" )
      ( "Points.PxConeR", po::value<VS>(), "PXCONE R value" )
      ( "Points.Donkersycutd", po::value<VS>(), "Donkers Ycut D values" )
      ( "Points.Donkersycutj", po::value<VS>(), "Donkers Ycut J values" )
      ( "Points.eminFractionPoints", po::value<VS>(), "Emin fraction points" )
      ( "Points.RPoints", po::value<VS>(), "R points" )
      ( "Points.eminFractionValue", po::value<VS>(), "Emin fraction value" )
      ( "Points.RValue", po::value<VS>(), "R value" );

    std::string filename= getItem<std::string>( "General.sjmGeneralOptions" );
    std::ifstream generalOptionsFile( filename );
    if( generalOptionsFile ) {
      po::store( po::parse_config_file( generalOptionsFile,
					generalOptions ), vm );
    }
    else {
      throw std::runtime_error( "General options file not found" );
    }
    setDefault( "General.url", "" );
    setDefault( "General.normalise", "1" );
    setDefault( "General.path", "." );
    setDefault( "Points.y34cut", "0.0045" );
    setDefault( "Points.y34y23cut", "0.5" );
    setDefault( "Points.PxConeEmin", "0.07" );
    setDefault( "Points.PxConeR", "0.7" );
    setDefault( "Points.eminFractionValue", "0.06" );
    setDefault( "Points.RValue", "0.7" );
    po::notify( vm );
    
  }
  catch( const po::error & ex ) {
    std::cerr << ex.what() << std::endl;
  }

  return;
}

// Set default values directly (builtin default(...) method
// of add_options(...) does not handle vector<string> argument)
void SjmConfigParser::setDefault( const std::string& key,
				  const std::string& value ) {
  if( not vm.count( key ) ) {
    vm.insert( std::make_pair( key,
			       po::variable_value( boost::any( VS { value } ), false ) ) );
  }
}

// Accessor template methods to options:

// Declare vector<string> specialisation so we can use it in
// generic template below:
template <>
VS SjmConfigParser::getItem( const std::string& tag ) const;

// For generic type convert lines via (special) operator>>:
template <typename T>
T SjmConfigParser::getItem( const std::string& tag ) const {
  VS lines= getItem<VS>( tag );
  T value;
  for( const std::string& line : lines ) {
    std::stringstream sstream( line );
    sstream >> value;
  }
  return value;
}

// For vector<string> return entry after checking its there,
// this is the basic retrieval of an option as vector<string>:
template <>
VS SjmConfigParser::getItem( const std::string& tag ) const {
  VS lines;
  if( vm.count( tag ) ) {
    lines= vm[tag].as<VS>();
  }
  else {
    std::ostringstream txt;
    txt << "SjmConfigParser::getItem: " << tag << " not found";
    throw std::runtime_error( txt.str() );
  }
  return lines;
}

// For single string return first entry of vector<string>:
template <>
std::string SjmConfigParser::getItem( const std::string& tag ) const {
  return getItem<VS>( tag )[0];
}

// Expect bool, int, float entries:
template bool SjmConfigParser::getItem( const std::string& tag ) const;
template int SjmConfigParser::getItem( const std::string& tag ) const;
template float SjmConfigParser::getItem( const std::string& tag ) const;

// For vector<double> et al define special operator>>:
template <typename T>
std::istream& operator>>( std::istream& input, std::vector<T>& vt ) {
  T t;
  while( input >> t ) vt.push_back( t );
  return input;
}
template std::vector<double> SjmConfigParser::getItem( const std::string& tag ) const;

// Filepath from path and file name:
VS SjmConfigParser::getFilepath( const std::string & tag ) const {
  VS files= getItem<VS>( tag );
  std::string path= getItem<std::string>( "General.path" );
  path+= "/";
  for( std::string & file : files ) {
    file.insert( 0, path );
  }
  return files;
}

// Points as vector<double> from "Points." section:
std::vector<double>
SjmConfigParser::getPoints( const std::string& tag ) const {
  std::string ptag= "Points."+tag;
  return getItem<std::vector<double>>( ptag );
}

// Print all configured items:
void SjmConfigParser::printConfig() const {
  std::vector< boost::shared_ptr<po::option_description> > cmdlineptions=
    cmdlineOptions.options();
  for( auto option : cmdlineptions ) {
    printOption( option );
  }
  std::vector< boost::shared_ptr<po::option_description> > specialopts=
    specialOptions.options();
  for( auto option : specialopts ) {
    printOption( option );
  }
  std::vector< boost::shared_ptr<po::option_description> > generalopts=
    generalOptions.options();
  for( auto option : generalopts ) {
    printOption( option );
  }
}

// Print a option:
void SjmConfigParser::printOption( boost::shared_ptr<po::option_description> option ) const {
  std::string key= option->key( "*" );
  cout << key << " (" << option->description() << "): ";
  if( vm.count( key ) ) {
    VS lines= getItem<VS>( key );
    if( lines.size() == 1 ) {
      cout << lines[0] << endl;
    }
    else {
      cout << endl;
      for( const std::string line : lines ) cout << line << endl;
    }
  }
  else {
    cout << "not given" << endl;
  }
  return;
}

