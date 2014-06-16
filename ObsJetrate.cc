
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

void ObsJetrate::getAndFillJetrateDataStructure( vector<Double_t> NJets,
						 Int_t Jetrate,
						 const string& tag ) {
  map<string,DataStructure*>::iterator iter= datastructures.find( tag );
  if( iter != datastructures.end() ) {
    JetrateDataStructure* jrds= 
      dynamic_cast<JetrateDataStructure*>( iter->second );
    if( jrds ) {
      jrds->fill( NJets, Jetrate );
    }
    else {
      std::cout << "ObsJetrate::getAndFillJetrateDataStructure: dynamic_cast to JetrateDataStructure failed for observable " 
		<< this->getName() << std::endl;
    }
  }
  else {
    std::cout << "ObsJetrate::getAndFillJetrateDataStructure: analysis " 
	      << tag << " not found" << std::endl;
  }
  return;
}

