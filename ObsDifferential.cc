
#include "ObsDifferential.hh"
#include "NtupleReader.hh"
#include "DataStructure.hh"
#include "DifferentialDataStructure.hh"

#include <iostream>
using std::cout;
using std::endl;

ObsDifferential::ObsDifferential( const string& name,
				  const vector<Double_t>& bins, 
				  const vector<Analysis>& variations ) :
  Observable(name), binedges(bins) {
  addAnalyses( variations, binedges );
}

void ObsDifferential::addAnalyses( const vector<Analysis>& variations, 
				   const vector<Double_t>& bins ) {
  for( size_t i= 0; i < variations.size(); i++ ) {      
    string tag= variations[i].getTag();
    datastructures[tag]= new DifferentialDataStructure( bins );
  }
}


void ObsDifferential::getAndFillDifferentialDataStructure( Double_t value, 
							   const string& tag ) {
  map<string,DataStructure*>::iterator iter= datastructures.find( tag );
  if( iter != datastructures.end() ) {
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
	      << tag << " not found" << std::endl;
  }
  return;
}
