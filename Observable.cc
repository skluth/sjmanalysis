
#include "Observable.hh"
#include "FilledObservable.hh"
#include "DataStructure.hh"
#include "MatrixDataStructure.hh"
#include <iostream>
using std::cout;
using std::endl;
#include<algorithm>

Observable::Observable( const string& namein ) : name( namein ) {}

Observable::~Observable() {}

void Observable::addAnalyses( const vector<Analysis>& va ) {
  analyses.insert( analyses.end(), va.begin(), va.end() );
  for( const Analysis& analysis: va ) {
    addAnalysis( analysis );
  }
}

class AnalysisEq {
  Analysis myanalysis;
public:
  AnalysisEq( const Analysis& analysis ) : myanalysis(analysis) {}
  bool operator()( const Analysis& analysis ) const {
    return myanalysis.getTag() == analysis.getTag();
  }
};
bool Observable::containsAnalysis( const Analysis& analysis ) {
  return std::find_if( analyses.begin(), analyses.end(), 
		       AnalysisEq(analysis) ) != analyses.end();
}

void Observable::fillAllAnalyses( NtupleReader* ntr ) {
  for( const Analysis& analysis: analyses ) {
    fill( ntr, analysis );
  }
}

void Observable::printVectorD( const string& txt, const vector<Double_t>& vd ) {
  cout << txt;
  for( size_t i= 0; i < vd.size(); i++ ) cout << " " << vd[i];
  cout << endl;
}

