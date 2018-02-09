
#include "ObservableFactory.hh"
#include "Observable.hh"

#include "ObsDifferential.hh"
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

ObservableFactory::ObservableFactory( const SjmConfigParser& sjmcp ) :
  sjmConfigs( sjmcp) {}

bool ObservableFactory::nameIs( const string& str, 
				const string& name ) {
  return str.find( name ) != string::npos;
}

// Handle all known observable names:
vector<Observable*>
ObservableFactory::createObservables( const vector<string>& obsnames,
				      const vector<Analysis>& analyses ) {
  vector<Observable*> vobs;
  for( const string & name : obsnames ) {
    if( nameIs( name, "thrust" ) ) {
      vobs.push_back( new ObsDifferential( "thrust", 
					   sjmConfigs.getPoints( "thrust" ),
					   analyses, 
					   new ThrustCalculator() ) );
    }
    else if( nameIs( name, "partonshower" ) ) {
      vobs.push_back( new ObsPartonShower( sjmConfigs.getPoints( "a14" ),
					   sjmConfigs.getPoints( "c202" ),
					   sjmConfigs.getPoints( "as" ),
					   sjmConfigs.getPoints( "mr" ),
					   analyses ) );
    }
    else if( nameIs( name, "durhamymerge23" ) ) {
      vobs.push_back( new ObsDifferential( "durhamymerge23", 
					   sjmConfigs.getPoints( "yNMPoints" ),
					   analyses,
					   new YnmdCalculator( 2 ) ) );
    }
    else if( nameIs( name, "jadeymerge23" ) ) {
      vobs.push_back( new ObsDifferential( "jadeymerge23", 
					   sjmConfigs.getPoints( "yNMPoints" ),
					   analyses,
					   new YnmjCalculator( 2 ) ) );
    }
    else if( nameIs( name, "durhamymergefj" ) ) {
      vobs.push_back( new ObsFastJetDiff( name, 
					  "eekt", 
					  sjmConfigs.getPoints( "yNMPoints" ),
					  analyses ) );
    }
    else if( nameIs( name, "jadeymergefj" ) ) {
      vobs.push_back( new ObsFastJetDiff( name, 
					  "jade", 
					  sjmConfigs.getPoints( "yNMPoints" ),
					  analyses ) );
    }
    else if( nameIs( name, "durhamycutfj" ) ) {
      vobs.push_back( new ObsJetrate( name, 
				      sjmConfigs.getPoints( "Donkersycutd" ),
				      analyses,
				      new FastJetYcutCalculator( "eekt" ) ) );
    }
    else if( nameIs( name, "jadeycutfj" ) ) {
      vobs.push_back( new ObsJetrate( name, 
				      sjmConfigs.getPoints( "Donkersycutj" ),
				      analyses,
				      new FastJetYcutCalculator( "jade" ) ) );
    }
    else if( nameIs( name, "durhamycut" ) ) {
      vobs.push_back( new ObsJetrate( name, 
				      sjmConfigs.getPoints( "Donkersycutd" ),
				      analyses,
				      new YcutCalculator( "durham" ) ) );
    }
    else if( nameIs( name, "jadeycut" ) ) {
      vobs.push_back( new ObsJetrate( name, 
				      sjmConfigs.getPoints( "Donkersycutj" ),
				      analyses,
				      new YcutCalculator( "jade" ) ) );
    }
    else if( nameIs( name, "antiktemin" ) ) {
      vobs.push_back( new ObsJetrate( name,
				      sjmConfigs.getPoints( "eminFractionPoints" ),
				      analyses,
				      new FastJetEminCalculator( "eeantikt", 
								 sjmConfigs.getItem<float>( "Points.RValue" ) ) ) );
    }
    else if( nameIs( name, "antiktR" ) ) {
      vobs.push_back( new ObsJetrate( name, 
				      sjmConfigs.getPoints( "RPoints" ),
				      analyses,
				      new FastJetRCalculator( "eeantikt", 
							      sjmConfigs.getItem<float>( "Points.eminFractionValue" ) ) ) );
    }
    else if( nameIs( name, "sisconeemin" ) ) {
      vobs.push_back( new ObsJetrate( name, 
				      sjmConfigs.getPoints( "eminFractionPoints" ),
				      analyses,
				      new FastJetEminCalculator( "eesiscone", 
								 sjmConfigs.getItem<float>( "Points.RValue" ) ) ) );
    }
    else if( nameIs( name, "sisconeR" ) ) {
      vobs.push_back( new ObsJetrate( name, 
				      sjmConfigs.getPoints( "RPoints" ),
				      analyses,
				      new FastJetRCalculator( "eesiscone", 
							      sjmConfigs.getItem<float>( "Points.eminFractionValue" ) ) ) );
    }
    else if( nameIs( name, "pxconeemin" ) ) 
      vobs.push_back( new ObsJetrate( name, 
				      sjmConfigs.getPoints( "PxEminPoints" ),
				      analyses,
				      new FastJetPxConeEminCalculator( sjmConfigs.getItem<float>( "Points.RValue" ) ) ) );
    else if( nameIs( name, "pxconeR" ) ) {
      vobs.push_back( new ObsJetrate( name, 
				      sjmConfigs.getPoints( "PxRPoints" ),
				      analyses,
				      new FastJetPxConeRCalculator( sjmConfigs.getItem<float>( "Points.PxConeEmin" ) ) ) );
    }
    else {
      string txt= "ObservableFactory::createObservables: wrong class name: " + name;
      throw std::logic_error( txt );
    }
  }
  cout << "ObservableFactory::createObservables:" << endl;
  for( size_t iobs= 0; iobs < vobs.size(); iobs++ ) {
    cout << vobs[iobs]->getName() << " ";
  }
  cout << endl;
  return vobs;
}

