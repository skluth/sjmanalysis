
#include "ObservableFactory.hh"
#include "Observable.hh"
#include "ObsThrust.hh"
#include "ObsPartonShower.hh"
#include "ObsDurhamYmerge23.hh"
#include "ObsJadeYmerge23.hh"
#include "ObsFastJetR.hh"
#include "ObsFastJetEmin.hh"
#include "ObsFastJetDiff.hh"
#include "ObsFastJetYcut.hh"
#include <iostream>
#include <algorithm>
using std::cout;
using std::endl;
#include <sstream>
using std::stringstream;
#include <cctype>

ObservableFactory::ObservableFactory() {

  // Thrust
  thrustbins.resize( 13 );
  thrustbins[0]= 0.00;
  thrustbins[1]= 0.01;
  thrustbins[2]= 0.02;
  thrustbins[3]= 0.03;
  thrustbins[4]= 0.04;
  thrustbins[5]= 0.05;
  thrustbins[6]= 0.07;
  thrustbins[7]= 0.09;
  thrustbins[8]= 0.12;
  thrustbins[9]= 0.15;
  thrustbins[10]= 0.22;
  thrustbins[11]= 0.30;
  thrustbins[12]= 0.50;

  // MR:
  mrbins.resize( 6 );
  mrbins[0]= 0.00;
  mrbins[1]= 0.06;
  mrbins[2]= 0.15;
  mrbins[3]= 0.38;
  mrbins[4]= 0.69;
  mrbins[5]= 1.00;

  // A14:
  for( size_t i= 0; i < 21; i++ ) a14bins.push_back( i*0.05 );

  // C202:
  c202bins.resize( 6 );
  c202bins[0]= 0.355;
  c202bins[1]= 0.402;
  c202bins[2]= 0.424;
  c202bins[3]= 0.446;
  c202bins[4]= 0.468;
  c202bins[5]= 0.495;

  // AS:
  asbins.resize( 8 );
  asbins[0]= -0.05;
  asbins[1]= 0.04;
  asbins[2]= 0.10;
  asbins[3]= 0.16;
  asbins[4]= 0.22;
  asbins[5]= 0.28;
  asbins[6]= 0.34;
  asbins[7]= 0.43;

  // incl. Jetrates and ynm distributions log scale:
  yNMbins.resize( 11 );
  for( size_t i= 0; i < yNMbins.size(); i++ ) {
    yNMbins[i]= 0.5*i;
  }

  // excl. jets Emin/Evis points:
  eminFraction.resize( 9 );
  eminFraction[0]= 0.02;
  eminFraction[1]= 0.04;
  eminFraction[2]= 0.06;
  eminFraction[3]= 0.08;
  eminFraction[4]= 0.10;
  eminFraction[5]= 0.12;
  eminFraction[6]= 0.14;
  eminFraction[7]= 0.16;
  eminFraction[8]= 0.18;

  // excl. jets R points:
  Rvalues.resize( 8 );
  Rvalues[0]= 0.2;
  Rvalues[1]= 0.4;
  Rvalues[2]= 0.6;
  Rvalues[3]= 0.7;
  Rvalues[4]= 0.8;
  Rvalues[5]= 1.0;
  Rvalues[6]= 1.2;
  Rvalues[7]= 1.4;

}

// Handle all known observable names:
vector<Observable*> ObservableFactory::createObservables( const vector<string>& obsnames,
							  const vector<Analysis>& analyses ) {
  vector<Observable*> vobs;
  for( size_t iobs= 0; iobs < obsnames.size(); iobs++ ) {
    string name= obsnames[iobs];
    if( name.find( "thrust" ) != string::npos )
      vobs.push_back( new ObsThrust( thrustbins, analyses ) );
    // else if( name == "mr" )
    //   vobs.push_back( new ObsMr( mrbins, analyses ) );
    else if( name.find( "partonshower" ) != string::npos )
      vobs.push_back( new ObsPartonShower( a14bins, 
					   c202bins,
					   asbins,
					   mrbins,
					   analyses ) );
    else if( name.find( "durhamymerge23" ) != string::npos ) 
      vobs.push_back( new ObsDurhamYmerge23( yNMbins, analyses ) );
    else if( name.find( "jadeymerge23" ) != string::npos ) 
      vobs.push_back( new ObsJadeYmerge23( yNMbins, analyses ) );
    else if( name.find( "durhamymergefj" ) != string::npos ) 
      vobs.push_back( new ObsFastJetDiff( name, "eekt", yNMbins, analyses ) );
    else if( name.find( "jadeymergefj" ) != string::npos )
      vobs.push_back( new ObsFastJetDiff( name, "jade", yNMbins, analyses ) );
    else if( name.find( "durhamycutfj" ) != string::npos ) 
      vobs.push_back( new ObsFastJetYcut( name, "eekt", yNMbins, analyses ) );
    else if( name.find( "jadeycutfj" ) != string::npos )
      vobs.push_back( new ObsFastJetYcut( name, "jade", yNMbins, analyses ) );
    else if( name.find( "antiktemin" ) != string::npos ) 
      vobs.push_back( new ObsFastJetEmin( name, "eeantikt", 0.7, eminFraction, analyses ) );
    else if( name.find( "antiktR" ) != string::npos ) 
      vobs.push_back( new ObsFastJetR( name, "eeantikt", 0.06, Rvalues, analyses ) );
    else if( name.find( "sisconeemin" ) != string::npos ) 
      vobs.push_back( new ObsFastJetEmin( name, "siscone", 0.7, eminFraction, analyses ) );
    else if( name.find( "sisconeR" ) != string::npos ) 
      vobs.push_back( new ObsFastJetR( name, "siscone", 0.06, Rvalues, analyses ) );
    else cout << "ObservableFactory::createObservables: name " << name
	      << " not recognised" << endl;
  }
  cout << "ObservableFactory::createObservables:" << endl;
  for( size_t iobs= 0; iobs < vobs.size(); iobs++ ) {
    cout << vobs[iobs]->getName() << " ";
  }
  cout << endl;
  return vobs;
}

