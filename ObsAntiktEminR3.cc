
#include "ObsAntiktEminR3.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"

ObsAntiktEminR3::ObsAntiktEminR3( const vector<Analysis>& variations, Double_t R ) : 
  Observable( "antikteminR3" ), Remin(R) {
  eminFraction.resize( 9 );
  eminFraction[0]= 0.02;
  eminFraction[1]= 0.04;
  eminFraction[2]= 0.06;
  eminFraction[3]= 0.08;
  eminFraction[4]= 0.1;
  eminFraction[5]= 0.12;
  eminFraction[6]= 0.14;
  eminFraction[7]= 0.16;
  eminFraction[8]= 0.18;
  addAnalyses( variations, eminFraction );
}

void ObsAntiktEminR3::fill( NtupleReader* ntr, const Analysis& variation ) {
  vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( variation.getReco() );
  size_t n= eminFraction.size();
  vector<Double_t> NJets( n );
  TFastJet tfjakt( vtlv, "eeantikt", Remin );
  Double_t Evis= tfjakt.Evis();
  for( size_t i= 0; i < n; i++ ) {
    vector<TLorentzVector> incljets= tfjakt.inclusive_jets( eminFraction[i]*Evis );
    NJets[i]= incljets.size();
  }
  datastructures[variation.getTag()].fill( NJets, 3 );
}
