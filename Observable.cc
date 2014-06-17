
#include "Observable.hh"
#include "FilledObservable.hh"
#include "DataStructure.hh"
#include <iostream>
using std::cout;
using std::endl;

Observable::Observable( const string& namein ) : name( namein ) {}

void Observable::print() const {
  for( map<string,DataStructure*>::const_iterator iter= datastructures.begin();
       iter != datastructures.end(); iter++ ) {
    cout << name << " " << iter->first << " events " << (iter->second)->getNEvents() << endl;
    (iter->second)->print();
  }
}

vector<FilledObservable*> Observable::getFilledObservables() const { 
  cout << "Observable::getFilledObservables: " << name 
       << ": create FilledObservable" << endl;
  FilledObservable* fobs= new FilledObservable( name, datastructures );
  vector<FilledObservable*> vfobs;
  vfobs.push_back( fobs );
  return vfobs;
}

