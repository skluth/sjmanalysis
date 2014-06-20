
#include "ObsDifferential.hh"
#include "NtupleReader.hh"
#include "DataStructure.hh"
#include "DifferentialDataStructure.hh"
#include "MatrixDataStructure.hh"

#include <iostream>
using std::cout;
using std::endl;

ObsDifferential::ObsDifferential( const string& name,
				  const vector<Double_t>& bins ) :
  Observable(name), binedges(bins) {}

// Call this method in derived class ctors:
void ObsDifferential::addAnalyses( const vector<Analysis>& variations ) {
  cout << "ObsDifferential::addAnalyses: adding analyses for " << getName() << endl;
  for( size_t ivar= 0; ivar < variations.size(); ivar++ ) {      
    string tag= variations[ivar].getTag();
    datastructures[tag]= new DifferentialDataStructure( binedges );
    if( variations[ivar].getReco2() != "none" ) {
      matrices[tag]= new MatrixDataStructure( binedges );
    }
  }
  return;
}

void 
ObsDifferential::getAndFillDifferentialDataStructure( Double_t value, 
						      const string& tag,
						      const map<string,DataStructure*>& dss,
						      Double_t weight ) {
  DataStructure* ds= getDataStructure( tag, dss );
  if( ds ) {
    DifferentialDataStructure* dds= getDifferentialDataStructure( ds );
    if( dds ) dds->fill( value, weight );
  }
  return;
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

