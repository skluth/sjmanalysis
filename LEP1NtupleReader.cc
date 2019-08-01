
#include "LEP1NtupleReader.hh"
#include <iostream>

LEP1NtupleReader::LEP1NtupleReader( const char* filename, 
				    const char* ntid, 
				    const bool lpr ) : 
  LEPNtupleReader( filename, ntid, lpr ) {
  SetBranchAddressChecked( "Itkmh", &nt_Itkmh );
}

bool LEP1NtupleReader::Preselection( const std::string& ecms ) {
  bool result= true;
  Int_t icjst= nt_Icjst;
  Int_t iebst= nt_Iebst;
  Int_t itkmh= nt_Itkmh;
  cutflow["cjst"]= icjst == 3;
  cutflow["ebst"]= iebst == 3;
  cutflow["tkmh"]= icjst == 3 and iebst == 3 and itkmh == 1;
  if( icjst != 3 or iebst != 3 or itkmh != 1 ) result= false;
  return result;
}

bool LEP1NtupleReader::Selection( const std::string& ecms ) {
  bool preselection= Preselection( ecms );
  bool nch7= nt_Ntkd02 >= 7;
  bool costt09= abscostt() < 0.9;
  cutflow["nch7"]= preselection and nch7;
  cutflow["costt09"]= preselection and costt09;
  bool result= preselection and nch7 and costt09;
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
