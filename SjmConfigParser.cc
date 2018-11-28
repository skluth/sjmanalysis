
#include "SjmConfigParser.hh"

#include <fstream>
#include <iostream>
#include <exception>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>
#include <boost/program_options.hpp>

namespace po= boost::program_options;

using std::cout;
using std::endl;

typedef std::vector<std::string> VS;

template <>
std::string SjmConfigParser::getItem( const std::string& tag ) const;

SjmConfigParser::SjmConfigParser( int argc, const char* argv[] ) {

  // Handle cmdline:
  try {
    po::options_description cmdlineOptions( "Cmdline" );
    cmdlineOptions.add_options()
      ( "help,h", "Help screen" )
      ( "config", po::value<VS>(), "Config file" );
    po::positional_options_description posOpts; 
    posOpts.add( "config", 1 );
    po::variables_map vm;
    po::store( po::command_line_parser( argc, argv ).
	       options( cmdlineOptions ).positional( posOpts ).run(),
	       vm );
    if( vm.count( "help" ) ) {
      std::cout << cmdlineOptions << std::endl;
    }
    if( not vm.count( "config" ) ) {
      throw std::runtime_error( "Config file name not given" );
    }    
    for( const auto & keyValue : vm ) {
      valuesMap[keyValue.first]= keyValue.second.as<VS>();
    }    
    descriptionMap["config"]= "Config file";
  }
  catch( const po::error & ex ) {
    std::cerr << ex.what() << std::endl;
  }

  // Handle config file from cmd line:
  std::map< std::string, std::string > specialOptsDescriptionMap;
  specialOptsDescriptionMap["General.test"]= "test field";
  specialOptsDescriptionMap["General.energy"]= "Energy point";
  specialOptsDescriptionMap["General.maxevt"]= "Max. number of events";
  specialOptsDescriptionMap["General.outfile"]= "Output file name";
  specialOptsDescriptionMap["General.sjmGeneralOptions"]= "Observable configuration file";
  specialOptsDescriptionMap["General.analyses"]= "Analyses list for LEP1, LEP15 or LEP2";
  specialOptsDescriptionMap["Data.lumi"]= "Data Lumi";
  specialOptsDescriptionMap["Data.files"]= "Data file names";
  specialOptsDescriptionMap["Signal.xsec"]= "Signal MC xsec";
  specialOptsDescriptionMap["Signal.name"]= "Signal MC name";
  specialOptsDescriptionMap["Signal.files"]= "Signal MC file names";
  specialOptsDescriptionMap["AltSignal.xsec"]= "Alt. signal MC xsec";
  specialOptsDescriptionMap["AltSignal.name"]= "Alt. signal MC name";
  specialOptsDescriptionMap["AltSignal.files"]= "Alt. signal MC file names";
  
  specialOptsDescriptionMap["BkgWWllqq.xsec"]= "Bkg WW->llqq MC xsec";
  specialOptsDescriptionMap["BkgWWllqq.name"]= "Bkg WW->llqq MC name";
  specialOptsDescriptionMap["BkgWWllqq.files"]= "Bkg WW->llqq MC file names";
  
  specialOptsDescriptionMap["BkgWWqqqq.xsec"]= "Bkg WW->qqqq MC xsec";
  specialOptsDescriptionMap["BkgWWqqqq.name"]= "Bkg WW->qqqq MC name";
  specialOptsDescriptionMap["BkgWWqqqq.files"]= "Bkg WW->qqqq MC file names";
  
  specialOptsDescriptionMap["BkgWWeeqq.xsec"]= "Bkg WW->eeqq MC xsec";
  specialOptsDescriptionMap["BkgWWeeqq.name"]= "Bkg WW->eeqq MC name";
  specialOptsDescriptionMap["BkgWWeeqq.files"]= "Bkg WW->eeqq MC file names";      
  
  valuesMap["General.test"]= VS { "bla" };
  // valuesMap["General.maxevt"]= VS { "999999" };
  valuesMap["General.maxevt"]= VS { "0" };
  valuesMap["General.sjmGeneralOptions"]= VS { "sjmGeneralOptions.cfg" };

  declareOptions( getItem<std::string>( "config" ),
		  specialOptsDescriptionMap );
    
  // General options:
  std::map< std::string, std::string > generalOptsDescriptionMap;
  generalOptsDescriptionMap["General.url"]= "Data URL";
  generalOptsDescriptionMap["General.normalise"]= "Normalise";
  generalOptsDescriptionMap["General.path"]= "Data path";
  
  generalOptsDescriptionMap["Observables.observable"]= "Observable names";
  generalOptsDescriptionMap["Observables.mtxunfold"]= "Mtx unfolding observables";
  
  generalOptsDescriptionMap["Points.thrust"]= "Thrust bin edges";      
  generalOptsDescriptionMap["Points.cpar"]= "C-parameter bin edges";      
  generalOptsDescriptionMap["Points.EEC"]= "EEC bin edges";
  
  generalOptsDescriptionMap["Points.y34cut"]= "y34 cut";
  generalOptsDescriptionMap["Points.y34y23cut"]= "y34/y23 cut";
  generalOptsDescriptionMap["Points.a14"]= "A14 bin edges";
  generalOptsDescriptionMap["Points.c202"]= "C202 bin edges";
  generalOptsDescriptionMap["Points.as"]= "AS bin edges";
  generalOptsDescriptionMap["Points.mr"]= "MR bin edges";
  
  generalOptsDescriptionMap["Points.yNMPoints"]= "Ynm points";
  generalOptsDescriptionMap["Points.PxEminPoints"]= "PXCONE Emin points";
  generalOptsDescriptionMap["Points.PxRPoints"]= "PXCONE R points";
  generalOptsDescriptionMap["Points.PxConeEmin"]= "PXCONE Emin value";
  generalOptsDescriptionMap["Points.PxConeR"]= "PXCONE R value";
  generalOptsDescriptionMap["Points.Donkersycutd"]= "Donkers Ycut D values";
  generalOptsDescriptionMap["Points.Donkersycutj"]= "Donkers Ycut J values";
  generalOptsDescriptionMap["Points.eminFractionPoints"]= "Emin fraction points";
  generalOptsDescriptionMap["Points.RPoints"]= "R points";
  generalOptsDescriptionMap["Points.eminFractionValue"]= "Emin fraction value";
  generalOptsDescriptionMap["Points.RValue"]= "R value";
  
  generalOptsDescriptionMap["LEP1Analyses.data"]= "LEP1 Data analyses";
  generalOptsDescriptionMap["LEP1Analyses.signal"]= "LEP1 Signal MC analyses";
  generalOptsDescriptionMap["LEP1Analyses.altsignal"]= "LEP1 Alt. signal MC analyses";

  generalOptsDescriptionMap["LEP15Analyses.data"]= "LEP15 Data analyses";
  generalOptsDescriptionMap["LEP15Analyses.signal"]= "LEP15 Signal MC analyses";
  generalOptsDescriptionMap["LEP15Analyses.altsignal"]= "LEP15 Alt. signal MC analyses";

  generalOptsDescriptionMap["LEP2Analyses.data"]= "LEP2 Data analyses";
  generalOptsDescriptionMap["LEP2Analyses.signal"]= "LEP2 Signal MC analyses";
  generalOptsDescriptionMap["LEP2Analyses.altsignal"]= "LEP2 Alt. signal MC analyses";
  generalOptsDescriptionMap["LEP2Analyses.bkgllqq"]= "LEP2 Bkg WW->llqq MC analyses";
  generalOptsDescriptionMap["LEP2Analyses.bkgqqqq"]= "LEP2 Bkg WW->qqqq MC analyses";
  generalOptsDescriptionMap["LEP2Analyses.bkgeeqq"]= "LEP2 Bkg WW->eeqq MC analyses";

  valuesMap["General.url"]= VS { "" };
  valuesMap["General.normalise"]= VS {"1" };
  valuesMap["General.path"]= VS { "." };
  valuesMap["Observables.mtxunfold"]= VS { "none" };
  valuesMap["Points.y34cut"]= VS { "0.0045" };
  valuesMap["Points.y34y23cut"]= VS { "0.5" };
  valuesMap["Points.PxConeEmin"]= VS { "0.07" };
  valuesMap["Points.PxConeR"]= VS { "0.7" };
  valuesMap["Points.eminFractionValue"]= VS { "0.06" };
  valuesMap["Points.RValue"]= VS { "0.7" };

  declareOptions( getItem<std::string>( "General.sjmGeneralOptions" ),
		  generalOptsDescriptionMap );

  return;
}

// Setup options described in map by reading file
void SjmConfigParser::
declareOptions( const std::string & filename,
		const std::map< std::string, std::string > & descriptions ) {
  try {
    po::options_description options( "Options" );
    for( const auto & keyValue : descriptions ) {
      options.add_options()( keyValue.first.c_str(),
			     po::value<VS>(),
			     keyValue.second.c_str() );
    }
    std::ifstream ifs( filename );
    po::variables_map vm;
    if( ifs ) po::store( po::parse_config_file( ifs, options ), vm );
    else throw std::runtime_error( "Config file not found" );
    po::notify( vm );    
    for( const auto & keyValue : vm ) {
      valuesMap[keyValue.first]= keyValue.second.as<VS>();
    }
    descriptionMap.insert( descriptions.begin(),
			   descriptions.end() );
  }
  catch( const po::error & ex ) {
    std::cerr << ex.what() << std::endl;
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
  if( valuesMap.count( tag ) ) {
    lines= valuesMap.at( tag );
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

// Print all options:
void SjmConfigParser::printConfig() const {
  cout << "SjmConfigParser: print all known options" << endl;
  for( auto keyValue : descriptionMap ) {
    cout << keyValue.first << " (" << keyValue.second << "): ";
    if( valuesMap.count( keyValue.first ) ) {
      VS lines= getItem<VS>( keyValue.first );
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
  }
}


