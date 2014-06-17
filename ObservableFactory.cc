
#include "ObservableFactory.hh"
#include "Observable.hh"
#include "ObsThrust.hh"
#include "ObsDurhamYmerge23.hh"
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
  yNMbins.resize( 11 );
  for( size_t i= 0; i < yNMbins.size(); i++ ) {
    yNMbins[i]= 0.5*i;
  }
}

// Handle all known observable names:
vector<Observable*> ObservableFactory::createObservables( const vector<string>& obsnames,
							  const vector<Analysis>& analyses ) {
  vector<Observable*> vobs;
  for( size_t iobs= 0; iobs < obsnames.size(); iobs++ ) {
    string name= obsnames[iobs];
    if( name == "thrust" )
      vobs.push_back( createThrust( analyses) );
    else if( name.find( "durhamymerge23" ) != string::npos ) 
      vobs.push_back( new ObsDurhamYmerge23( yNMbins, analyses ) );
    else if( name.find( "durhamymergefj" ) != string::npos ) 
      vobs.push_back( new ObsFastJetDiff( name, "eekt", yNMbins, analyses ) );
    else if( name.find( "jadeymergefj" ) != string::npos )
      vobs.push_back( new ObsFastJetDiff( name, "jade", yNMbins, analyses ) );
    else if( name.find( "durhamycutfj" ) != string::npos ) 
      vobs.push_back( createFastJetYcut( name, "eekt", analyses) );
    else if( name.find( "jadeycutfj" ) != string::npos )
      vobs.push_back( createFastJetYcut( name, "jade", analyses) );
    else if( name.find( "antiktemin" ) != string::npos ) 
      vobs.push_back( createFastJetEmin( name, "eeantikt", analyses ) );
    else if( name.find( "antiktR" ) != string::npos ) 
      vobs.push_back( createFastJetR( name, "eeantikt", analyses ) );
    else if( name.find( "sisconeemin" ) != string::npos ) 
      vobs.push_back( createFastJetEmin( name, "siscone", analyses ) );
    else if( name.find( "sisconeR" ) != string::npos ) 
      vobs.push_back( createFastJetR( name, "siscone", analyses ) );
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

// Thrust:
Observable* ObservableFactory::createThrust( const vector<Analysis>& analyses ) {
  vector<Double_t> thrustbins( 13 );
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
  ObsThrust* obsthrust= new ObsThrust( thrustbins, analyses );
  return obsthrust;
}

// Jetrates at fixed y_cut points, expect e.g. 3 at end of name to signal 3-jet rate:
Observable* ObservableFactory::createFastJetYcut( const string& obsname, 
						  const string& algo,
						  const vector<Analysis>& analyses ) {
  vector<Double_t> ycutpoints( 11 );
  for( size_t i= 0; i < ycutpoints.size(); i++ ) {
    ycutpoints[i]= 0.5*i;
  }
  //Int_t njet= getNjetFromName( obsname );
  return new ObsFastJetYcut( obsname, algo, ycutpoints, analyses );
}



// N-jet rate vs Emin/Evis for fixed R, expect e.g. 3 at end of name to signal
// 3-jet rate:
Observable* ObservableFactory::createFastJetEmin( const string& obsname,
						  const string& algo,
						  const vector<Analysis>& analyses,
						  Double_t rvalue ) {
  vector<Double_t> eminFraction( 9 );
  eminFraction[0]= 0.02;
  eminFraction[1]= 0.04;
  eminFraction[2]= 0.06;
  eminFraction[3]= 0.08;
  eminFraction[4]= 0.10;
  eminFraction[5]= 0.12;
  eminFraction[6]= 0.14;
  eminFraction[7]= 0.16;
  eminFraction[8]= 0.18;
  //Int_t njet= getNjetFromName( obsname );
  ObsFastJetEmin* obsfastjetemin= new ObsFastJetEmin( obsname, algo, rvalue, 
						      eminFraction, analyses );
  return obsfastjetemin;
}

// N-jet rate vs R for fixed Emin/Evis, convention as above:
Observable* ObservableFactory::createFastJetR( const string& obsname, 
					       const string& algo,
					       const vector<Analysis>& analyses,
					       Double_t eminfrac ) {
  vector<Double_t> Rvalues( 8 );
  Rvalues[0]= 0.2;
  Rvalues[1]= 0.4;
  Rvalues[2]= 0.6;
  Rvalues[3]= 0.7;
  Rvalues[4]= 0.8;
  Rvalues[5]= 1.0;
  Rvalues[6]= 1.2;
  Rvalues[7]= 1.4;
  //  Int_t njet= getNjetFromName( obsname );
  //  if( njet > 0 and njet < 7 ) {
  ObsFastJetR* obsfastjetR= new ObsFastJetR( obsname, algo, eminfrac, 
					     Rvalues, analyses );
    //}
  return obsfastjetR;
}

// Extract number from last field of string name:
Int_t ObservableFactory::getNjetFromName( const string& name ) {
  Int_t njet= 0;
  size_t intpos= name.size()-1;
  if( isdigit( name[intpos] ) ) {
    stringstream sname( name.substr( intpos, 1 ) );
    sname >> njet;
  }
  else {
    cout << "ObservableFactory::getNjetFromName: not an integer in last field " 
	 << name[intpos] << ", " << name << endl;      
  }
  return njet;
}

