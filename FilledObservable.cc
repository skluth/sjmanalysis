
#include "FilledObservable.hh"

#include "DataStructure.hh"
#include "DifferentialDataStructure.hh"
#include "JetrateDataStructure.hh"
#include "MatrixDataStructure.hh"

#include <iostream>
using std::cout;
using std::endl;

using std::map;
using std::string;

FilledObservable::FilledObservable( const string& obsname,
				    const map<string,DifferentialDataStructure*>& ddss,
				    const map<string,MatrixDataStructure*>& mdss ) : 
  name(obsname), lerrorMatrices(false), migrationMatrices(mdss) {
  for( const auto & keyValue : ddss ) datastructures[keyValue.first]= keyValue.second;
}

FilledObservable::FilledObservable( const string& obsname,
				    const map<string,JetrateDataStructure*>& jrdss ) :
  name(obsname), lerrorMatrices(false) {
  for( const auto & keyValue : jrdss ) datastructures[keyValue.first]= keyValue.second;
}

void FilledObservable::finalise() {
  for( map<string,DataStructure*>::iterator iter= datastructures.begin();
       iter != datastructures.end(); iter++ ) {
    (iter->second)->normalise();
  }
}

void FilledObservable::Print() const {
  // for( map<string,DataStructure*>::const_iterator iter= datastructures.begin();
  //      iter != datastructures.end(); iter++ ) {
  //   cout << name << " " << iter->first << " events " << (iter->second)->getNEvents() << endl;
  //   (iter->second)->print();
  // }
  for( const auto & keyValue : datastructures ) {
    cout << name << " " << keyValue.first << " events " << (keyValue.second)->getNEvents()
	 << endl;
    (keyValue.second)->Print();
  }
  if( migrationMatrices.empty() ) {
    cout << "No migration matrices" << endl;
  }
  else {
    for( map<string,MatrixDataStructure*>::const_iterator iter= migrationMatrices.begin();
	 iter != migrationMatrices.end(); iter++ ) {
      cout << name << " " << iter->first << " events " << (iter->second)->getNEvents() << endl;
      (iter->second)->Print();
    }
  }
}

const map<string,DataStructure*>& FilledObservable::getData() const { 
  return datastructures;
}
const map<string,MatrixDataStructure*>& FilledObservable::getMigrationMatrices() const { 
  return migrationMatrices;
}
map<string,MatrixDataStructure*> FilledObservable::getErrorMatrices() const { 
  map<string,MatrixDataStructure*> errorMatrices;
  for( const auto & keyValue : datastructures ) {
    if( keyValue.second != 0 ) {
      errorMatrices[keyValue.first]= keyValue.second->getErrorMatrix();
    }
  }
  return errorMatrices;
}

DataStructure* FilledObservable::getDataStructure( const Analysis& anal ) const {
  string tag= anal.getTag();
  map<string,DataStructure*>::const_iterator iter= datastructures.find( tag );
  DataStructure* ds= 0;
  if( iter != datastructures.end() ) ds= iter->second;
  else cout << "FilledObservable::getDataStructure: " << tag << " not found" << endl;
  return ds;
}

void FilledObservable::setDataStructure( DataStructure* dsp, const Analysis& anal ) {
  string tag= anal.getTag();
  datastructures[tag]= dsp;
}


bool FilledObservable::containsAnalysis( const Analysis& anal ) const {
  string tag= anal.getTag();
  bool result= true;
  if( datastructures.find( tag ) == datastructures.end() ) result= false;
  return result;
}

