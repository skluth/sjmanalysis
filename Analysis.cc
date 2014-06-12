
#include "Analysis.hh"

#include <iostream>
using std::cout;
using std::endl;

Analysis::Analysis( const string& source, const string& reco, const string& cuts, 
		    const string& mccuts,  const string& reco2,
		    const string& unfoldsource, const string& unfoldmethod ) :
  a_source(source), a_reco(reco), a_cuts(cuts), a_mccuts(mccuts), a_reco2(reco2),
  a_unfoldsource(unfoldsource), a_unfoldmethod(unfoldmethod) {}

Analysis::Analysis() {}
Analysis::~Analysis() {}

string Analysis::getSource() const { return a_source; }
string Analysis::getReco() const { return a_reco; }
string Analysis::getCuts() const { return a_cuts; }
string Analysis::getMccuts() const { return a_mccuts; }
string Analysis::getReco2() const { return a_reco2; }
string Analysis::getUnfoldSource() const { return a_unfoldsource; }
string Analysis::getUnfoldMethod() const { return a_unfoldmethod; }

string Analysis::getTag() const { 
  string result= a_source + " " + a_reco + " " + a_cuts;
  if( a_mccuts != "none" ) {
    result+= " ";
    result+= a_mccuts;
  }
  if( a_reco2 != "none" ) {
    result+= " ";
    result+= a_reco2;
  }
  if( a_unfoldsource != "none" ) {
    result+= " ";
    result+= a_unfoldsource;
  }
  if( a_unfoldmethod != "none" ) {
    result+= " ";
    result+= a_unfoldmethod;
  }
  return result;
}

void Analysis::print() const {
  cout << a_source << " " << a_reco << " " << a_cuts << " " << a_mccuts << " " 
       << a_reco2 << " " << a_unfoldsource << " " << a_unfoldmethod << endl;
}

