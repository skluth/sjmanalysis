
#include "HepMC2AnalysisProcessor.hh"

#include <iostream>
#include <stdexcept>

HepMC2AnalysisProcessor::HepMC2AnalysisProcessor( const std::string & filename ) :
  reader( filename ) {}

void HepMC2AnalysisProcessor::runAnalysis() {

  for( int ievent= 0; ievent < 1000; ievent++ ) {

    if( reader.GetEvent() ) {
      
      const std::vector<TLorentzVector> vtlvisr= reader.GetLorentzVectors( "isr" );

      if( reader.MCNonRad() ) {
      
      const std::vector<TLorentzVector> vtlvh= reader.GetLorentzVectors( "hadron" );
      std::cout << "HepMC2AnalysisProcessor::runAnalysis: event " << ievent+1 << std::endl;
      std::cout << "Hadron level" << std::endl;
      TLorentzVector sumh;
      for( const TLorentzVector & tlv : vtlvh ) {
	std::cout << tlv.Px() << " "
		  << tlv.Py() << " "
		  << tlv.Pz() << " "
		  << tlv.E() << " "
		  << tlv.M() << std::endl;
	sumh+= tlv;
      }
      std::cout << "4-vector sum hadron level" << std::endl;
      std::cout << sumh.Px() << " "
		  << sumh.Py() << " "
		  << sumh.Pz() << " "
		  << sumh.E() << " "
		  << sumh.M() << std::endl;

      for( const TLorentzVector & tlv : vtlvisr ) {
	sumh+= tlv;
      }
      std::cout << "4-vector sum hadron level incl isr" << std::endl;
      std::cout << sumh.Px() << " "
		  << sumh.Py() << " "
		  << sumh.Pz() << " "
		  << sumh.E() << " "
		  << sumh.M() << std::endl;
            
      const std::vector<TLorentzVector> vltvp= reader.GetLorentzVectors( "parton" );
      std::cout << "Parton level" << std::endl;
      TLorentzVector sump;
      for( const TLorentzVector & tlv : vltvp ) {
	std::cout << tlv.Px() << " "
		  << tlv.Py() << " "
		  << tlv.Pz() << " "
		  << tlv.E() << " "
		  << tlv.M() << std::endl;
	sump+= tlv;
      }
      std::cout << "4-vector sum parton level" << std::endl;
      std::cout << sump.Px() << " "
		  << sump.Py() << " "
		  << sump.Pz() << " "
		  << sump.E() << " "
		  << sump.M() << std::endl;
      
      for( const TLorentzVector & tlv : vtlvisr ) {
	sump+= tlv;
      }
      std::cout << "4-vector sum parton level incl isr" << std::endl;
      std::cout << sumh.Px() << " "
		  << sumh.Py() << " "
		  << sumh.Pz() << " "
		  << sumh.E() << " "
		  << sumh.M() << std::endl;
      }
      else {
	std::cout << "HepMC2AnalysisProcessor: MC radiative event " << ievent+1
		  << std::endl;
	for( const TLorentzVector & tlv : vtlvisr ) tlv.Print();
      }
             
    }
    else {
      std::cout << "HepMC2AnalysisProcessor::runAnalysis: read " << ievent
		<< " events" << std::endl;
      break;
      //      throw std::runtime_error( "HepMC2Reader::GetEvent() failure" );
    }
  }
    
  
}
