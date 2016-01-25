
#include "FilledObservable.hh"
#include "DataStructure.hh"
#include "DifferentialDataStructure.hh"
#include "MatrixDataStructure.hh"
#include <iostream>
using std::cout;
using std::endl;

FilledObservable::FilledObservable( const string& obsname,
				    const map<string,DataStructure*>& dss,
				    const map<string,MatrixDataStructure*>& mds ) : 
  name(obsname), datastructures(dss), matrices(mds) {}

FilledObservable::FilledObservable( const string& obsname,
				    const map<string,DifferentialDataStructure*>& ddss,
				    const map<string,MatrixDataStructure*>& mds ) : 
  name(obsname), matrices(mds) {
  for( map<string,DifferentialDataStructure*>::const_iterator iter= ddss.begin();
       iter != ddss.end(); iter++ ) {
    datastructures[iter->first]= iter->second;
  }
}

void FilledObservable::finalise() {
  for( map<string,DataStructure*>::iterator iter= datastructures.begin();
       iter != datastructures.end(); iter++ ) {
    (iter->second)->normalise();
  }
}

void FilledObservable::print() const {
  for( map<string,DataStructure*>::const_iterator iter= datastructures.begin();
       iter != datastructures.end(); iter++ ) {
    cout << name << " " << iter->first << " events " << (iter->second)->getNEvents() << endl;
    (iter->second)->print();
  }
  if( matrices.empty() ) {
    cout << "No matrices" << endl;
  }
  else {
    for( map<string,MatrixDataStructure*>::const_iterator iter= matrices.begin();
	 iter != matrices.end(); iter++ ) {
      cout << name << " " << iter->first << " events " << (iter->second)->getNEvents() << endl;
      (iter->second)->print();
    }
  }
}

const map<string,DataStructure*>& FilledObservable::getData() const { 
  return datastructures;
}
const map<string,MatrixDataStructure*>& FilledObservable::getMatrices() const { 
  return matrices;
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

bool FilledObservable::containsAnalysis( const Analysis& anal ) {
  string tag= anal.getTag();
  bool result= true;
  if( datastructures.find( tag ) == datastructures.end() ) result= false;
  return result;
}

