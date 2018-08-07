
#include "YcutCalculator.hh"

#include "NtupleReader.hh"
#include "TMath.h"
#include <iostream>

using std::string;
using std::vector;

YcutCalculator::YcutCalculator( const string & algo ) : 
  algorithm(algo) {}

void YcutCalculator::print() const {
  std::cout << "YcutCalculator algorithm: " << algorithm << std::endl;
}

vector<Double_t> 
YcutCalculator::getValues( NtupleReader* ntr, 
			   const vector<Double_t> & Ycutpoints,
			   const string & reco ) const {
  size_t n= Ycutpoints.size();
  vector<Double_t> NJets( n );
  for( size_t i= 0; i < n; i++ ) {
    // Double_t ycut= TMath::Power( 10.0, -Ycutpoints[i] );
    Double_t ycut= Ycutpoints[i];
    for( int njet= 2; njet <= 6; njet++ ) {
      Double_t ymergelo, ymergehi;
      if( algorithm == "jade" ) {
	ymergehi= ntr->getYmergeE( reco, njet-1 );
	ymergelo= ntr->getYmergeE( reco, njet );
      }
      else if( algorithm == "durham" ) {
	ymergehi= ntr->getYmergeD( reco, njet-1 );
	ymergelo= ntr->getYmergeD( reco, njet );
      }
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

