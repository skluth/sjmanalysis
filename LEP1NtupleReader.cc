
#include "LEP1NtupleReader.hh"
#include <iostream>

LEP1NtupleReader::LEP1NtupleReader( const char* filename, 
				    const char* ntid, 
				    const bool lpr ) : 
  NtupleReader( filename, ntid, lpr ) {
  SetBranchAddressChecked( "Itkmh", &nt_Itkmh );
}

bool LEP1NtupleReader::Preselection( const std::string& ecms ) {
  bool result= true;
  Int_t icjst= nt_Icjst;
  Int_t iebst= nt_Iebst;
  Int_t itkmh= nt_Itkmh;
  if( icjst != 3 or iebst != 3 or itkmh != 1 ) result= false;
  return result;
}

bool LEP1NtupleReader::Selection( const std::string& ecms ) {
  bool result= true;
  if( ! Preselection( ecms ) ) result= false;
  // if( nt_Ntkd02 < 5 ) result= false;
  if( nt_Ntkd02 < 7 ) result= false;
  if( abscostt() > 0.9 ) result= false;
  return result;
}

const std::map<std::string,bool> 
LEP1NtupleReader::getSelections( const std::string& ecms ) {
  std::map<std::string,bool> selections;
  selections["stand"]= false;
  selections["costt07"]= false;
  // selections["nch7"]= false;
  bool standardSelection= Selection( ecms );
  if( standardSelection ) {
    selections["stand"]= true;
    // if( nt_Ntkd02 >= 7 ) selections["nch7"]= true;
    if( abscostt() < 0.7 ) selections["costt07"]= true;
  }
  return selections;
}
