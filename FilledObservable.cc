
#include "FilledObservable.hh"
#include "DataStructure.hh"
#include <iostream>
using std::cout;
using std::endl;

FilledObservable::FilledObservable( const string& namein,
				    const map<string,DataStructure*>& dssin ) : 
  name(namein), datastructures(dssin) {}

void FilledObservable::finalise() {
  for( map<string,DataStructure*>::iterator iter= datastructures.begin();
       iter != datastructures.end(); iter++ ) {
    (iter->second)->normalise();
  }
}

void FilledObservable::print() {
  for( map<string,DataStructure*>::iterator iter= datastructures.begin();
       iter != datastructures.end(); iter++ ) {
    cout << name << " " << iter->first << " events " << (iter->second)->getNEvents() << endl;
    (iter->second)->print();
  }
}

//string FilledObservable::getName() { return name; }

map<string,DataStructure*> FilledObservable::getData() { 
  map<string,DataStructure*> result( datastructures.begin(), datastructures.end() );
  return result; 
}

// void FilledObservable::setData( const map<string,DataStructure*>& dssin) {
//   datastructures= dssin;
// }


DataStructure* FilledObservable::getDataStructure( const Analysis& anal ) {
  string tag= anal.getTag();
  return datastructures[tag];
}

void FilledObservable::setDataStructure( DataStructure* dsp, const Analysis& anal ) {
  string tag= anal.getTag();
  datastructures[tag]= dsp;
}
