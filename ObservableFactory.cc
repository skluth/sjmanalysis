
#include "ObservableFactory.hh"

#include "Observable.hh"
#include "ObsDifferential.hh"
#include "ObsEEC.hh"
#include "ThrustCalculator.hh"
#include "YnmdCalculator.hh"
#include "YnmjCalculator.hh"

#include "ObsPartonShower.hh"
#include "ObsFastJetDiff.hh"

#include "ObsJetrate.hh"
#include "YcutCalculator.hh"
#include "FastJetYcutCalculator.hh"
#include "FastJetEminCalculator.hh"
#include "FastJetRCalculator.hh"
#include "FastJetPxConeEminCalculator.hh"
#include "FastJetPxConeRCalculator.hh"

#include "SjmConfigParser.hh"

#include <iostream>
#include <algorithm>
using std::cout;
using std::endl;
#include <sstream>
using std::stringstream;
#include <cctype>
#include <stdexcept>

using std::string;
using std::vector;

ObservableFactory::ObservableFactory( const SjmConfigParser & sjmcp ) :
  sjmConfigs( sjmcp) {}

bool ObservableFactory::nameIs( const string & str, 
				const string & name ) {
  // return str.find( name ) != string::npos;
  return str == name;
}

// Handle all known observable names:
vector<Observable*>
ObservableFactory::createObservables( const vector<string> & obsnames,
				      const vector<Analysis> & analyses ) {
  vector<Observable*> vobs;
  vector<Double_t> ynmpoints= sjmConfigs.getPoints( "yNMPoints" );
  // first (0.0) and last ynm are (should be(!)) outside of kinematic range,
  // so don't run jet rates there
  vector<Double_t> ycutpoints( ynmpoints.begin()+1, ynmpoints.end()-2 );
  for( const string & name : obsnames ) {
    Observable* obsp= 0;
    if( nameIs( name, "thrust" ) ) {
      obsp= new ObsDifferential( "thrust", 
				 sjmConfigs.getPoints( "thrust" ),
				 analyses, 
				 new ThrustCalculator() );
    }
    else if( nameIs( name, "EECnsc" ) ) {
      obsp= new ObsEEC( name,
			sjmConfigs.getPoints( "EEC" ),
			analyses,
			false );
    }
    else if( nameIs( name, "EEC" ) ) {
      obsp= new ObsEEC( name,
			sjmConfigs.getPoints( "EEC" ),
			analyses );
    }
    else if( nameIs( name, "partonshower" ) ) {
      obsp= new ObsPartonShower( sjmConfigs.getPoints( "a14" ),
				 sjmConfigs.getPoints( "c202" ),
				 sjmConfigs.getPoints( "as" ),
				 sjmConfigs.getPoints( "mr" ),
				 analyses );
    }
    else if( nameIs( name, "durhamymerge23" ) ) {
      obsp= new ObsDifferential( name, 
				 ynmpoints,
				 analyses,
				 new YnmdCalculator( 2 ) );
    }
    else if( nameIs( name, "jadeymerge23" ) ) {
      obsp= new ObsDifferential( name,
				 ynmpoints,
				 analyses,
				 new YnmjCalculator( 2 ) );
    }
    else if( nameIs( name, "durhamymergefj" ) ) {
      obsp= new ObsFastJetDiff( name, 
				"eekt", 
				ynmpoints,
				analyses );
    }
    else if( nameIs( name, "jadeymergefj" ) ) {
      obsp= new ObsFastJetDiff( name, 
				"jade", 
				ynmpoints,
				analyses );
    }
    else if( nameIs( name, "durhamycutfj" ) ) {
      obsp= new ObsJetrate( name, 
			    // sjmConfigs.getPoints( "Donkersycutd" ),
			    ycutpoints,
			    analyses,
			    new FastJetYcutCalculator( "eekt" ) );
    }
    else if( nameIs( name, "jadeycutfj" ) ) {
      obsp= new ObsJetrate( name, 
			    // sjmConfigs.getPoints( "Donkersycutj" ),
			    ycutpoints,
			    analyses,
			    new FastJetYcutCalculator( "jade" ) );
    }
    else if( nameIs( name, "durhamycut" ) ) {
      obsp= new ObsJetrate( name, 
			    // sjmConfigs.getPoints( "Donkersycutd" ),
			    ycutpoints,
			    analyses,
			    new YcutCalculator( "durham" ) );
    }
    else if( nameIs( name, "jadeycut" ) ) {
      obsp= new ObsJetrate( name, 
			    // sjmConfigs.getPoints( "Donkersycutj" ),
			    ycutpoints,
			    analyses,
			    new YcutCalculator( "jade" ) );
    }
    else if( nameIs( name, "antiktemin" ) ) {
      obsp= new ObsJetrate( name,
			    sjmConfigs.getPoints( "eminFractionPoints" ),
			    analyses,
			    new FastJetEminCalculator( "eeantikt", 
						       sjmConfigs.getItem<float>( "Points.RValue" ) ) );
    }
    else if( nameIs( name, "antiktR" ) ) {
      obsp= new ObsJetrate( name, 
			    sjmConfigs.getPoints( "RPoints" ),
			    analyses,
			    new FastJetRCalculator( "eeantikt", 
						    sjmConfigs.getItem<float>( "Points.eminFractionValue" ) ) );
    }
    else if( nameIs( name, "sisconeemin" ) ) {
      obsp= new ObsJetrate( name, 
			    sjmConfigs.getPoints( "eminFractionPoints" ),
			    analyses,
			    new FastJetEminCalculator( "eesiscone", 
						       sjmConfigs.getItem<float>( "Points.RValue" ) ) );
    }
    else if( nameIs( name, "sisconeR" ) ) {
      obsp= new ObsJetrate( name, 
			    sjmConfigs.getPoints( "RPoints" ),
			    analyses,
			    new FastJetRCalculator( "eesiscone", 
						    sjmConfigs.getItem<float>( "Points.eminFractionValue" ) ) );
    }
    else if( nameIs( name, "pxconeemin2" ) ) {
      obsp= new ObsJetrate( name, 
			    sjmConfigs.getPoints( "eminFractionPoints" ),
			    analyses,
			    new FastJetEminCalculator( "pxcone",
						       sjmConfigs.getItem<float>( "Points.RValue" ) ) );
    }
    else if( nameIs( name, "pxconeR2" ) ) {
      obsp= new ObsJetrate( name,
			    sjmConfigs.getPoints( "RPoints" ),
			    analyses,
			    new FastJetRCalculator( "pxcone",
						    sjmConfigs.getItem<float>( "Points.eminFractionValue" ) ) );
    }
    else if( nameIs( name, "pxconeemin" ) ) {
      obsp= new ObsJetrate( name, 
			    sjmConfigs.getPoints( "PxEminPoints" ),
			    analyses,
			    new FastJetPxConeEminCalculator( sjmConfigs.getItem<float>( "Points.PxConeR" ) ) );
    }
    else if( nameIs( name, "pxconeR" ) ) {
      obsp= new ObsJetrate( name, 
			    sjmConfigs.getPoints( "PxRPoints" ),
			    analyses,
			    new FastJetPxConeRCalculator( sjmConfigs.getItem<float>( "Points.PxConeEmin" ) ) );
    }    
    else {
      string txt= "ObservableFactory::createObservables: wrong class name: " + name;
      throw std::logic_error( txt );
    }
    vobs.push_back( obsp );
  }
  cout << "ObservableFactory::createObservables:" << endl;
  for( const Observable* obsp : vobs ) cout << obsp->getName() << " ";
  cout << endl;
  return vobs;
}

