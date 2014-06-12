
#include "ObsJetrate.hh"
#include "DataStructure.hh"
#include "JetrateDataStructure.hh"
#include <iostream>
using std::cout;
using std::endl;

ObsJetrate::ObsJetrate( string name ) : Observable( name ) {}

void ObsJetrate::addAnalyses( const vector<Analysis>& variations, 
			      const vector<Double_t>& points ) {
  for( size_t i= 0; i < variations.size(); i++ ) {      
    string tag= variations[i].getTag();
    datastructures[tag]= new JetrateDataStructure( points );
  }
}

