
#include "Observable.hh"
#include "FilledObservable.hh"
#include "DataStructure.hh"
#include "MatrixDataStructure.hh"
#include <iostream>
using std::cout;
using std::endl;

Observable::Observable( const string& namein ) : name( namein ) {}

Observable::~Observable() {
  deleteDataStructures( datastructures );
}

void Observable::deleteDataStructures( map<string,DataStructure*>& dss ) {
  for( map<string,DataStructure*>::iterator iter= dss.begin();
       iter != dss.end(); iter++ ) {
    DataStructure* ds= iter->second;
    if( ds ) delete ds;
  }
}

void Observable::print() const {
  for( map<string,DataStructure*>::const_iterator iter= datastructures.begin();
       iter != datastructures.end(); iter++ ) {
    cout << name << " " << iter->first << " events " << (iter->second)->getNEvents() << endl;
    (iter->second)->print();
  }
  for( map<string,MatrixDataStructure*>::const_iterator iter= matrices.begin();
       iter != matrices.end(); iter++ ) {
    cout << name << " " << iter->first << " events " << (iter->second)->getNEvents() << endl;
    (iter->second)->print();
  }
}

vector<FilledObservable*> Observable::getFilledObservables() const { 
  cout << "Observable::getFilledObservables: " << name 
       << ": create FilledObservable" << endl;
  FilledObservable* fobs= new FilledObservable( name, datastructures, matrices );
  vector<FilledObservable*> vfobs;
  vfobs.push_back( fobs );
  return vfobs;
}

DataStructure* 
Observable::getDataStructure( const string& tag,
			      const map<string,DataStructure*>& dss ) const {
  DataStructure* ds= 0;
  map<string,DataStructure*>::const_iterator iter= dss.find( tag );
  if( iter != dss.end() ) {
    ds= iter->second;
  }
  else {
    cout << "Observable::getDataStructure: analysis " 
	 << tag << " not found for " << name << endl;
  }
  return ds;
}

bool Observable::containsAnalysis( const Analysis& anal ) {
  return containsAnalysisInDataStructure( anal, datastructures );
}

bool Observable::containsAnalysisInDataStructure( const Analysis& anal,
						  const map<string,DataStructure*>& dss ) {
  string tag= anal.getTag();
  bool result= true;
  if( dss.find( tag ) == dss.end() ) result= false;
  return result;
}

void Observable::printVectorD( const string& txt, const vector<Double_t>& vd ) {
  cout << txt;
  for( size_t i= 0; i < vd.size(); i++ ) cout << " " << vd[i];
  cout << endl;
}

