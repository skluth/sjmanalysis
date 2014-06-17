
#include "ObsJetrate.hh"
#include "DataStructure.hh"
#include "JetrateDataStructure.hh"
#include "FilledObservable.hh"
#include <iostream>
using std::cout;
using std::endl;

ObsJetrate::ObsJetrate( string name, const vector<Double_t>& pts ) : 
  Observable( name ), points(pts) {}

void ObsJetrate::addAnalyses( const vector<Analysis>& variations ) {
  for( size_t i= 0; i < variations.size(); i++ ) {      
    string tag= variations[i].getTag();
    jetrates2[tag]= new JetrateDataStructure( points, 2 );
    jetrates3[tag]= new JetrateDataStructure( points, 3 );
    jetrates4[tag]= new JetrateDataStructure( points, 4 );
    jetrates5[tag]= new JetrateDataStructure( points, 5 );
    jetrates6[tag]= new JetrateDataStructure( points, 6 );
  }
}

void ObsJetrate::getAndFillJetrateDataStructures( const vector<Double_t>& NJets,
						  const string& tag ) {
  getAndFillJetrateDataStructure( NJets, tag, jetrates2 );
  getAndFillJetrateDataStructure( NJets, tag, jetrates3 );
  getAndFillJetrateDataStructure( NJets, tag, jetrates4 );
  getAndFillJetrateDataStructure( NJets, tag, jetrates5 );
  getAndFillJetrateDataStructure( NJets, tag, jetrates6 );
  return;
}

void ObsJetrate::getAndFillJetrateDataStructure( const vector<Double_t>& NJets,
						 const string& tag,
						 map<string,DataStructure*>& jetrates ) {
  map<string,DataStructure*>::iterator iter= jetrates.find( tag );
  if( iter != jetrates.end() ) {
    JetrateDataStructure* jrds= 
      dynamic_cast<JetrateDataStructure*>( iter->second );
    if( jrds ) {
      jrds->fill( NJets );
    }
    else {
      std::cout << "ObsJetrate::getAndFillJetrateDataStructure: dynamic_cast to JetrateDataStructure failed for observable " 
		<< this->getName() << std::endl;
    }
  }
  else {
    std::cout << "ObsJetrate::getAndFillJetrateDataStructures: analysis " 
	      << tag << " not found" << std::endl;
  }
  return; 
}

vector<FilledObservable*> ObsJetrate::getFilledObservables() const {
  cout << "ObsJetrate::getFilledObservables: " << name << ": create FilledObservables" << endl;
  FilledObservable* fobs2jetrate= new FilledObservable( name+"R2", jetrates2 );
  FilledObservable* fobs3jetrate= new FilledObservable( name+"R3", jetrates3 );
  FilledObservable* fobs4jetrate= new FilledObservable( name+"R4", jetrates4 );
  FilledObservable* fobs5jetrate= new FilledObservable( name+"R5", jetrates5 );
  FilledObservable* fobs6jetrate= new FilledObservable( name+"R6", jetrates6 );
  vector<FilledObservable*> vfobs;
  vfobs.push_back( fobs2jetrate );
  vfobs.push_back( fobs3jetrate );
  vfobs.push_back( fobs4jetrate );
  vfobs.push_back( fobs5jetrate );
  vfobs.push_back( fobs6jetrate );
  return vfobs;
}

void ObsJetrate::print() const {
  printDatastructures( jetrates2 );
  printDatastructures( jetrates3 );
  printDatastructures( jetrates4 );
  printDatastructures( jetrates5 );
  printDatastructures( jetrates6 );
}

void ObsJetrate::printDatastructures( const map<string,DataStructure*>& dss ) const {  
  for( map<string,DataStructure*>::const_iterator iter= dss.begin();
       iter != dss.end(); iter++ ) {
    cout << name << " " << iter->first << " events " << (iter->second)->getNEvents() << endl;
    (iter->second)->print();
  }
}

