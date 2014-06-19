
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
				const vector<Double_t>& bins, 
				const vector<Analysis>& variations ) :
  ObsDifferential( name, bins ), Algorithm( algo ) {
  addAnalyses( variations );
}

ObsFastJetDiff::~ObsFastJetDiff() {
  deleteDataStructures( ymerge23 );
  deleteDataStructures( ymerge34 );
  deleteDataStructures( ymerge45 );
  deleteDataStructures( ymerge56 );
}

void ObsFastJetDiff::addAnalyses( const vector<Analysis>& variations ) {
  cout << "ObsFastJetDiff::addAnalyses: adding analyses for " << getName() << endl;
  for( size_t i= 0; i < variations.size(); i++ ) {      
    string tag= variations[i].getTag();
    ymerge23[tag]= new DifferentialDataStructure( binedges );
    ymerge34[tag]= new DifferentialDataStructure( binedges );
    ymerge45[tag]= new DifferentialDataStructure( binedges );
    ymerge56[tag]= new DifferentialDataStructure( binedges );
  }
}

void ObsFastJetDiff::fill( NtupleReader* ntr, const Analysis& variation ) {
  const vector<TLorentzVector>& vtlv= ntr->GetLorentzVectors( variation.getReco() );
  TFastJet tfj( vtlv, Algorithm.c_str() );
  Double_t yflip= tfj.ymerge( 2 );
  getAndFillDifferentialDataStructure( -TMath::Log10( yflip ), variation.getTag(),
				       ymerge23 );
  yflip= tfj.ymerge( 3 );
  getAndFillDifferentialDataStructure( -TMath::Log10( yflip ), variation.getTag(),
				       ymerge34 );
  yflip= tfj.ymerge( 4 );
  getAndFillDifferentialDataStructure( -TMath::Log10( yflip ), variation.getTag(),
				       ymerge45 );
  yflip= tfj.ymerge( 5 );
  getAndFillDifferentialDataStructure( -TMath::Log10( yflip ), variation.getTag(),
				       ymerge56 );
  return;
}

vector<FilledObservable*> ObsFastJetDiff::getFilledObservables() const {
  cout << "ObsFastJetDiff::getFilledObservables: " << name 
       << ": create FilledObservables" << endl;  
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

bool ObsFastJetDiff::containsAnalysis( const Analysis& anal ) {
  return containsAnalysisInDataStructure( anal, ymerge23 );
}

