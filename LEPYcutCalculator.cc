
#include "LEPYcutCalculator.hh"

#include "NtupleReader.hh"
#include "LEPNtupleReader.hh"
#include "TMath.h"
#include <iostream>
#include <stdexcept>

using std::string;
using std::vector;

LEPYcutCalculator::LEPYcutCalculator( const string & algo ) : 
  algorithm(algo) {}

void LEPYcutCalculator::print() const {
  std::cout << "LEPYcutCalculator algorithm: " << algorithm << std::endl;
}

vector<Double_t> 
LEPYcutCalculator::getValues( NtupleReader* ntr, 
			      const vector<Double_t> & Ycutpoints,
			      const string & reco ) const {
  LEPNtupleReader* lepntr= dynamic_cast<LEPNtupleReader*>( ntr );
  if( lepntr == nullptr ) {
    throw std::runtime_error( "LEPYcutCalculator::getValues: no LEP ntuple" );
  }
  size_t n= Ycutpoints.size();
  vector<Double_t> NJets( n );
  for( size_t i= 0; i < n; i++ ) {
    Double_t ycut= Ycutpoints[i];
    for( int njet= 2; njet <= 6; njet++ ) {
      Double_t ymergehi= lepntr->getYmerge( algorithm, reco, njet-1 );
      Double_t ymergelo= lepntr->getYmerge( algorithm, reco, njet );
      NJets[i]= 0.0;
      if( ymergelo > 0.0 and ymergehi > 0.0 and 
	  ymergelo < ycut and ycut < ymergehi ) {
	NJets[i]= njet;
	break;
      }
    }
  }
  return NJets;
}

