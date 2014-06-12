
#include "ObsAntiktRR3.hh"
#include "NtupleReader.hh"
#include "TFastJet.hh"

ObsAntiktRR3::ObsAntiktRR3( const vector<Analysis>& variations, Double_t eminfrac ) : 
  Observable::Observable( "antiktRR3" ), EminFraction(eminfrac) {
  Rvalues.resize( 8 );
  Rvalues[0]= 0.2;
  Rvalues[1]= 0.4;
  Rvalues[2]= 0.6;
  Rvalues[3]= 0.7;
  Rvalues[4]= 0.8;
  Rvalues[5]= 1.0;
  Rvalues[6]= 1.2;
  Rvalues[7]= 1.4;
  addAnalyses( variations, Rvalues );
}

void ObsAntiktRR3::fill( NtupleReader* ntr, const Analysis& variation ) {
  vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( variation.getReco() );
  size_t n= Rvalues.size();
  vector<Double_t> NJets( n );
  for( size_t i= 0; i < n; i++ ) {
    TFastJet tfjakt( vtlv, "eeantikt", Rvalues[i] );
    Double_t Evis= tfjakt.Evis();
    vector<TLorentzVector> incljets= tfjakt.inclusive_jets( EminFraction*Evis );
    NJets[i]= incljets.size();
  }
  datastructures[variation.getTag()].fill( NJets, 3 );
}

