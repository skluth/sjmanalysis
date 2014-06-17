
#include "ObsDifferential.hh"
#include "NtupleReader.hh"
#include "DataStructure.hh"
#include "DifferentialDataStructure.hh"

#include <iostream>
using std::cout;
using std::endl;

ObsDifferential::ObsDifferential( const string& name,
				  const vector<Double_t>& bins ) :
  Observable(name), binedges(bins) {
  //  addAnalyses( variations, binedges );
}

void ObsDifferential::addAnalyses( const vector<Analysis>& variations ) {
  cout << "ObsDifferential::addAnalyses: adding analyses for " << getName() << endl;
  for( size_t i= 0; i < variations.size(); i++ ) {      
    string tag= variations[i].getTag();
    datastructures[tag]= new DifferentialDataStructure( binedges );
  }
}

void ObsDifferential::getAndFillDifferentialDataStructure( Double_t value, 
							   const string& tag,
							   map<string,DataStructure*>& dss ) {
  map<string,DataStructure*>::iterator iter= dss.find( tag );
  if( iter != dss.end() ) {
    DifferentialDataStructure* dds= 
      dynamic_cast<DifferentialDataStructure*>( iter->second );
    if( dds ) {
      dds->fill( value );
    }
    else {
      std::cout << "ObsDifferential::getAndFillDifferentialDataStructure: dynamic_cast to DifferentialDataStructure failed for observable " 
		<< this->getName() << std::endl;
    }
  }
  else {
    std::cout << "ObsDifferential::getAndFillDifferentialDataStructure: analysis " 
	      << tag << " not found for " << this->getName() << std::endl;
  }
  return;
}
