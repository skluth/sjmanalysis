
#include "Observable.hh"
#include "DataStructure.hh"
#include <iostream>
using std::cout;
using std::endl;

Observable::Observable( string namein ) : name( namein ) {}

void Observable::finalise() {
  map<string,DataStructure*> data= getData();
  for( map<string,DataStructure*>::iterator iter= data.begin();
       iter != data.end(); iter++ ) {
    (iter->second)->normalise();
  }
}

void Observable::print() {
  map<string,DataStructure*> data= getData();
  for( map<string,DataStructure*>::iterator iter= data.begin();
       iter != data.end(); iter++ ) {
    cout << name << " " << iter->first << " events " << (iter->second)->getNEvents() << endl;
    (iter->second)->print();
  }
}

string Observable::getName() { return name; }

map<string,DataStructure*> Observable::getData() { 
  map<string,DataStructure*> result( datastructures.begin(), datastructures.end() );
  return result; 
}

DataStructure* Observable::getDataStructure( const Analysis& anal ) {
  string tag= anal.getTag();
  return datastructures[tag];
}

void Observable::setDataStructure( DataStructure* dsp, const Analysis& anal ) {
  string tag= anal.getTag();
  datastructures[tag]= dsp;
}
