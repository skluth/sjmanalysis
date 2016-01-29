
#include "ObsFastJetDiff.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"
#include "DataStructure.hh"
#include "DifferentialDataStructure.hh"
#include "FilledObservable.hh"

#include <iostream>
using std::cout;
using std::endl;
#include "TMath.h"

ObsFastJetDiff::ObsFastJetDiff( const string& name, const string& algo,
				const vector<Double_t>& ynmbins, 
				const vector<Analysis>& variations,
				const bool lprint ) :
  Observable( name ), Algorithm( algo ), binedges(ynmbins) {
  addAnalyses( variations );
  if( lprint ) {
    cout << "ObsFastJetDiff::ObsFastJetDiff: create " << getName() 
	 << " with algorithm " << algo << " for y23, y34, y45, y56" << endl;
    printVectorD( "Binedges:", ynmbins );
  }
}

ObsFastJetDiff::~ObsFastJetDiff() {}

//void ObsFastJetDiff::addAnalyses( const vector<Analysis>& variations ) {
void ObsFastJetDiff::addAnalysis( const Analysis& analysis ) {
  //  for( size_t i= 0; i < variations.size(); i++ ) {      
  // string tag= variations[i].getTag();
  string tag= analysis.getTag();
  ymerge23[tag]= new DifferentialDataStructure( binedges );
  ymerge34[tag]= new DifferentialDataStructure( binedges );
  ymerge45[tag]= new DifferentialDataStructure( binedges );
  ymerge56[tag]= new DifferentialDataStructure( binedges );
  //  }
}

void ObsFastJetDiff::fill( NtupleReader* ntr, const Analysis& variation ) {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( variation.getReco() );
  TFastJet tfj( vtlv, Algorithm.c_str() );
  string tag= variation.getTag();
  ymerge23.at(tag)->fill( -TMath::Log10( tfj.ymerge( 2 ) ) );
  ymerge34.at(tag)->fill( -TMath::Log10( tfj.ymerge( 3 ) ) );
  ymerge45.at(tag)->fill( -TMath::Log10( tfj.ymerge( 4 ) ) );
  ymerge56.at(tag)->fill( -TMath::Log10( tfj.ymerge( 5 ) ) );
  return;
}

vector<FilledObservable*> ObsFastJetDiff::getFilledObservables() const {
  //cout << "ObsFastJetDiff::getFilledObservables: " << name 
  //     << ": create FilledObservables" << endl;  
  FilledObservable* fobsymerge23= new FilledObservable( name+"23", ymerge23 );
  FilledObservable* fobsymerge34= new FilledObservable( name+"34", ymerge34 );
  FilledObservable* fobsymerge45= new FilledObservable( name+"45", ymerge45 );
  FilledObservable* fobsymerge56= new FilledObservable( name+"56", ymerge56 );
  vector<FilledObservable*> vfobs;
  vfobs.push_back( fobsymerge23 );
  vfobs.push_back( fobsymerge34 );
  vfobs.push_back( fobsymerge45 );
  vfobs.push_back( fobsymerge56 );
  return vfobs;
}


