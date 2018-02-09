
#include "Analysis.hh"

#include <iostream>
using std::cout;
using std::endl;
#include <boost/algorithm/string.hpp>
#include <exception>

Analysis::Analysis( const string& src, const string& r, const string& c, 
		    const string& mcc,  const string& r2,
		    const string& bkgsts,
		    const string& unfsrc, const string& unfm ) :
  source(src), reco(r), cuts(c), mccuts(mcc), reco2(r2),
  bkgstatus(bkgsts),
  unfoldsource(unfsrc), unfoldmethod(unfm) {}

Analysis::Analysis( const string& options ) : 
  source(), reco(), cuts(), mccuts( "none" ), reco2( "none" ),
  bkgstatus( "none" ),
  unfoldsource( "none" ), unfoldmethod( "none" ) {
  std::vector<std::string> tokens;
  boost::split( tokens, options, boost::is_any_of( " " ), boost::token_compress_on );
  size_t ntoken= tokens.size();
  if( ntoken < 3 ) throw std::runtime_error( "Analysis::Analysis: options wrong" );
  source= tokens[0];
  reco= tokens[1];
  cuts= tokens[2];
  if( ntoken > 3 ) mccuts= tokens[3];
  if( ntoken > 4 ) reco2= tokens[4];
  if( ntoken > 5 ) bkgstatus= tokens[5];
  if( ntoken > 6 ) unfoldsource= tokens[6];
  if( ntoken > 7 ) unfoldmethod= tokens[7];
}


Analysis::Analysis() {}

Analysis::~Analysis() {}

string Analysis::getSource() const { return source; }
string Analysis::getReco() const { return reco; }
string Analysis::getCuts() const { return cuts; }
string Analysis::getMccuts() const { return mccuts; }
string Analysis::getReco2() const { return reco2; }
string Analysis::getBkgStatus() const { return bkgstatus; }
string Analysis::getUnfoldSource() const { return unfoldsource; }
string Analysis::getUnfoldMethod() const { return unfoldmethod; }

void Analysis::setSource( const string& src ) { source=src; }
void Analysis::setReco( const string& r ) { reco=r; }
void Analysis::setCuts( const string& c ) { cuts=c; }
void Analysis::setMccuts( const string& mcc ) { mccuts=mcc; }
void Analysis::setReco2( const string& r2 ) { reco2=r2; }
void Analysis::setBkgStatus( const string& bkgsts ) { bkgstatus=bkgsts; }
void Analysis::setUnfoldSource( const string& unfsrc ) { unfoldsource=unfsrc; }
void Analysis::setUnfoldMethod( const string& unfm ) { unfoldmethod=unfm; }

string Analysis::getTag() const { 
  string result= source + " " + reco + " " + cuts;
  if( mccuts != "none" ) {
    result+= " ";
    result+= mccuts;
  }
  if( reco2 != "none" ) {
    result+= " ";
    result+= reco2;
  }
  if( bkgstatus != "none" ) {
    result+= " ";
    result+= bkgstatus;
  }
  if( unfoldsource != "none" ) {
    result+= " ";
    result+= unfoldsource;
  }
  if( unfoldmethod != "none" ) {
    result+= " ";
    result+= unfoldmethod;
  }
  return result;
}

void Analysis::print() const {
  cout << source << " " << reco << " " << cuts << " " << mccuts << " " 
       << reco2 << " " << bkgstatus << " "
       << unfoldsource << " " << unfoldmethod << endl;
}

