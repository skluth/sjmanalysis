
#include "ObsGroomed.hh"

#include "NtupleReader.hh"
#include "DifferentialDataStructure.hh"
#include "MatrixDataStructure.hh"
#include "FilledObservable.hh"
#include "TLorentzVector.h"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"

#include <iostream>
#include <cmath>

// Andriis opalJGR
// #include "Groom.h"


ObsGroomed::ObsGroomed( const std::vector<Double_t> & grthbins,
			const std::vector<Double_t> & grcpbins,
			const std::vector<Analysis> & variations,
			bool lprint ) :
  Observable( "groomedshapes" ),
  binedges { { "grthrust", grthbins }, { "grcpar", grcpbins } },
  betaValues { 0.0, 1.0 },
  zcutValues { 0.05, 0.10, 0.15 } {
  addAnalyses( variations );
  if( lprint ) {
    std::cout << "ObsGroomed::ObsGroomed: create " << getName() << std::endl;
    printVectorD( "Groomed thrust binedges:", binedges.at( "grthrust" ) );
    printVectorD( "Groomed C-parameter binedges:", binedges.at( "grcpar" ) );
    printVectorD( "Grooming beta values:", betaValues );
    printVectorD( "Grooming zcut values:", zcutValues );
  }
}

void ObsGroomed::addAnalysis( const Analysis & analysis ) {
  std::string tag= analysis.getTag();
  std::string reco2= analysis.getReco2();
  LoopFunc func= [this,&tag,&reco2] ( const std::string & key,
				      const std::string & bzkey ) {
    DdsMap & ddsmap= grData[bzkey];
    ddsmap[tag]= new DifferentialDataStructure( binedges.at( key ) );
    MdsMap & mdsmap= grMatrices[bzkey];
    if( reco2 != "none" ) {
      mdsmap[tag]= new MatrixDataStructure( binedges.at( key ) );
    }    
  };
  loop( func );
  return;
}

std::string makeKey( const std::string & prefix,
		     const Double_t beta, const Double_t zcut ) {
  std::string strbeta= std::to_string( beta );
  std::string strzcut= std::to_string( zcut );
  return prefix + "_" + strbeta.substr( 0, 3 ) + "_" + strzcut.substr( 0, 4 );
}
void ObsGroomed::loop( const LoopFunc & func ) const {
  for( Double_t beta : betaValues ) {
    for( Double_t zcut : zcutValues ) {
      for( const std::string & key : { "grthrust", "grcpar" } ) {
	std::string bzkey= makeKey( key, beta, zcut );
	func( key, bzkey );
      }
    }
  }
  return;
}

Double_t dot3( const fastjet::PseudoJet & pj1,
	       const fastjet::PseudoJet & pj2 ) {
  return pj1.px()*pj2.px() + pj1.py()*pj2.py() + pj1.pz()*pj2.pz();
}
Double_t jetThrust( const fastjet::PseudoJet & jet,
		    const std::vector<fastjet::PseudoJet> & constituents ) {
  Double_t jetT= 0.0;
  for( const fastjet::PseudoJet & constituent : constituents ) {
    jetT+= std::abs( dot3( jet, constituent ) );
  }
  jetT/= jet.modp();
  return jetT;
}
const std::map<std::string,Double_t>
getValue( const std::vector<fastjet::PseudoJet> & cajets,
	  const Double_t beta, const Double_t zcut ) {

  // Setup softdrop jet transformer:
  Double_t R0= 1.0;
  fastjet::contrib::SoftDrop softdrop( beta, zcut,
				       fastjet::contrib::RecursiveSymmetryCutBase::theta_E,
				       R0,
				       std::numeric_limits<Double_t>::infinity(),
				       fastjet::contrib::RecursiveSymmetryCutBase::larger_E );

  // Apply softdrop to the two C/A jets and get their constituents:
  fastjet::PseudoJet cajet1= cajets[0];
  fastjet::PseudoJet cajet2= cajets[1];
  fastjet::PseudoJet sdjet1= softdrop( cajet1 );
  fastjet::PseudoJet sdjet2= softdrop( cajet2 );
  std::vector<fastjet::PseudoJet> sdconstituents1= sdjet1.constituents();
  std::vector<fastjet::PseudoJet> sdconstituents2= sdjet2.constituents();
  std::vector<fastjet::PseudoJet> sdconstituents;
  sdconstituents.insert( sdconstituents.end(), sdconstituents1.begin(), sdconstituents1.end() );
  sdconstituents.insert( sdconstituents.end(), sdconstituents2.begin(), sdconstituents2.end() );  
  
  // Groomed thrust:
  Double_t T1= jetThrust( cajet1, sdconstituents1 );
  Double_t T2= jetThrust( cajet2, sdconstituents2 );
  Double_t sump= 0.0;
  for( const fastjet::PseudoJet & sdconstituent : sdconstituents ) sump+= sdconstituent.modp();
  Double_t grThrust= 1.0 - ( T1 + T2 ) / sump;

  // Groomed C-parameter:
  Double_t grCpar= 0.0;
  for( const fastjet::PseudoJet & sdconstituent1 : sdconstituents ) {
    Double_t p1= sdconstituent1.modp();
    for( const fastjet::PseudoJet & sdconstituent2 : sdconstituents ) {
      Double_t p2= sdconstituent2.modp();
      Double_t p1p2= p1*p2;
      Double_t sdc1sdc2op1p2= dot3( sdconstituent1, sdconstituent2 )/p1p2;
      grCpar+= p1p2*( 1.0 - sdc1sdc2op1p2*sdc1sdc2op1p2 );
    }
  }
  grCpar= 3.0/2.0*grCpar / (sump*sump);
  
  // Return results:
  std::map<std::string,Double_t> result;
  result[ makeKey( "grthrust", beta, zcut ) ]= grThrust;
  result[ makeKey( "grcpar", beta, zcut ) ]= grCpar;
  return result;

}

const std::map<std::string,Double_t>
ObsGroomed::getValues( NtupleReader* ntr,
		       const std::string & reco ) const {
  std::vector<TLorentzVector> vtlv= ntr->GetLorentzVectors( reco );
  std::map<std::string,Double_t> result;
  if( vtlv.size() >= 2 ) {
    fastjet::JetDefinition jet_def( fastjet::ee_genkt_algorithm, 3.0*3.14159/2.0, 0.0 );
    fastjet::ClusterSequence clusseq( vtlv, jet_def );
    std::vector<fastjet::PseudoJet> cajets= clusseq.exclusive_jets( 2 );
    for( Double_t beta : betaValues ) {
      for( Double_t zcut : zcutValues ) {
	const std::map<std::string,Double_t> value= getValue( cajets, beta, zcut );
	result.insert( value.begin(), value.end() );

	// from Andriis opalJGR
	// std::vector<fastjet::PseudoJet> CPARTONS1;
	// std::vector<fastjet::PseudoJet> CPARTONS2;
	// fastjet::PseudoJet AXPARTONS1;
	// fastjet::PseudoJet AXPARTONS2;
	// GroomInput<TLorentzVector>( vtlv, beta, zcut,
	// 			    CPARTONS1, CPARTONS2,
	// 			    AXPARTONS1, AXPARTONS2 );
	// std::map<std::string,double>
	//   obsPARTONS= GetGroomedObservables<fastjet::PseudoJet>( CPARTONS1,
	// 							 CPARTONS2,
	// 							 AXPARTONS1,
	// 							 AXPARTONS2 );
	// std::string grthkey= makeKey( "grthrust", beta, zcut );
	// std::string grcpkey= makeKey( "grcpar", beta, zcut );
	// std::cout << beta << " " << zcut << " "
	// 	  << result.at( grthkey ) << " " << 1.0-obsPARTONS["T"] << " "
	// 	  << result.at( grcpkey )  << " " << obsPARTONS["C"]
	// 	  << std::endl;
	
      }
    }
    
  }
  else {
    std::cout << "ObsGroomed::getValues: less than two objects, return -1" << std::endl;
    LoopFunc func= [this,&result] ( const std::string & key, const std::string & bzkey ) {
      result[bzkey]= -1.0;
    };
    loop( func );
  }
  return result;
}

void ObsGroomed::fill( NtupleReader* ntr, const Analysis & variation ) { 
  std::string reco= variation.getReco();
  std::string tag= variation.getTag();
  const std::map<std::string,Double_t> values= getValues( ntr, reco );
  LoopFunc func= [this,&tag,&values] ( const std::string & key, const std::string & bzkey ) {
    const DdsMap & ddsmap= grData.at( bzkey );
    ddsmap.at( tag )->fill( values.at( bzkey ) );
  };
  loop( func );
  std::string reco2= variation.getReco2();
  if( reco2 != "none" and ntr->isMC() ) {
    const std::map<std::string,Double_t> MCvalues= getValues( ntr, reco2 );
    LoopFunc func2= [this,&tag,&MCvalues,&values] ( const std::string & key,
						    const std::string & bzkey ) {
      const MdsMap & mdsmap= grMatrices.at( bzkey );
      mdsmap.at( tag )->fill( MCvalues.at( bzkey ), values.at( bzkey ) );
    };
    loop( func2 );
  }
  return;
}

std::vector<FilledObservable*> ObsGroomed::getFilledObservables() const {
  std::cout << "ObsGroomed::getFilledObservables: " 
	    << "create groomed thrust" << std::endl;
  std::vector<FilledObservable*> vfobs;
  LoopFunc func= [this,&vfobs] ( const std::string & key, const std::string & bzkey ) {
    const DdsMap & ddsmap= grData.at( bzkey );      
    const MdsMap & mdsmap= grMatrices.at( bzkey );
    vfobs.push_back( new FilledObservable( bzkey, ddsmap, mdsmap ) );
  };
  loop( func );
  return vfobs;
}
