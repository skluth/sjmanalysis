
#include "ObsJetrate.hh"
#include "DataStructure.hh"
#include "JetrateDataStructure.hh"
#include "FilledObservable.hh"
#include "JetrateCalculator.hh"

#include <iostream>
using std::cout;
using std::endl;

ObsJetrate::ObsJetrate( string name, 
			const vector<Double_t>& pts,
			const vector<Analysis>& variations, 
			const JetrateCalculator* calc,
			const bool lprint ) :
  Observable( name ), points( pts ), calculator( calc ) {
  addAnalyses( variations );
  if( lprint ) {
    cout << "ObsJetrate::ObsJetrate: 2, 3, 4, 5, 6-jet fractions for " << name << endl;
    printVectorD( "Points:", pts );
  }
}

ObsJetrate::~ObsJetrate() {}

//void ObsJetrate::addAnalyses( const vector<Analysis>& variations ) {
void ObsJetrate::addAnalysis( const Analysis& analysis ) {
  //  for( size_t i= 0; i < variations.size(); i++ ) {
  //    string tag= variations[i].getTag();
  string tag= analysis.getTag();
  jetrates2[tag]= new JetrateDataStructure( points, 2 );
  jetrates3[tag]= new JetrateDataStructure( points, 3 );
  jetrates4[tag]= new JetrateDataStructure( points, 4 );
  jetrates5[tag]= new JetrateDataStructure( points, 5 );
  jetrates6[tag]= new JetrateDataStructure( points, 6 );
    //  }
}

void ObsJetrate::fill( NtupleReader* ntr, const Analysis& variation ) {
  vector<Double_t> NJets= calculator->getValues( ntr,
						 points,
						 variation.getReco() );
  string tag= variation.getTag();
  jetrates2.at(tag)->fill( NJets );
  jetrates3.at(tag)->fill( NJets );
  jetrates4.at(tag)->fill( NJets );
  jetrates5.at(tag)->fill( NJets );
  jetrates6.at(tag)->fill( NJets );
  return;
}

vector<FilledObservable*> ObsJetrate::getFilledObservables() const {
  // cout << "ObsJetrate::getFilledObservables: " << name << ": create FilledObservables" << endl;
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

void ObsJetrate::printDatastructures( const map<string,JetrateDataStructure*>& jrds ) const {  
  for( map<string,JetrateDataStructure*>::const_iterator iter= jrds.begin();
       iter != jrds.end(); iter++ ) {
    cout << name << " " << iter->first << " events " << (iter->second)->getNEvents() << endl;
    (iter->second)->print();
  }
}

