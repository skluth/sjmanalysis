
#include "ObsDifferential.hh"
#include "NtupleReader.hh"
#include "DataStructure.hh"
#include "DifferentialDataStructure.hh"

#include <iostream>
using std::cout;
using std::endl;

ObsDifferential::ObsDifferential( const string& name,
				  const vector<Double_t>& bins ) :
  Observable(name), binedges(bins) {}

// Call in derived class ctors
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
  DataStructure* ds= getDataStructure( tag, dss );
  if( ds ) {
    DifferentialDataStructure* dds= getDifferentialDataStructure( ds );
    if( dds ) dds->fill( value );
  }
}

DifferentialDataStructure* 
ObsDifferential::getDifferentialDataStructure( DataStructure* ds ) const {
  DifferentialDataStructure* dds= dynamic_cast<DifferentialDataStructure*>( ds );
  if( not dds ) {
    cout << "ObsDifferential::getDifferentialDataStructure: dynamic_cast to DifferentialDataStructure failed for observable " 
	 << name << endl;
  }
  return dds;
}

