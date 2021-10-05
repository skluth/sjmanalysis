
#include "NtupleReader.hh"
#include <algorithm>

#include <numeric>

#include <iostream>
#include <numeric>

const std::map<std::string,bool> NtupleReader::getCutflow() {
  return cutflow;
}

Double_t Evis( const std::vector<TLorentzVector> & v ) {
  Double_t evis= std::accumulate( v.begin(), v.end(), 0.0,
                                  []( Double_t sum, const TLorentzVector & tlv ) {
                                    return sum+= tlv.E();
                                  }
                                  );
  return evis;
}

