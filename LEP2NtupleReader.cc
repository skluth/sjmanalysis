
#include "LEP2NtupleReader.hh"
#include <iostream>
#include <map>

LEP2NtupleReader::LEP2NtupleReader( const char* filename, 
				    const char* ntid, 
				    const bool lpr ) : 
  LEPNtupleReader( filename, ntid, lpr ),
  ecmsranges{
    { "130", range( 129.0, 131.0 ) },
    { "136", range( 135.0, 137.0 ) },
    { "161", range( 160.0, 162.0 ) },
    { "172", range( 170.0, 173.0 ) },
    { "183", range( 180.0, 184.5 ) },
    { "189", range( 188.0, 190.0 ) },
    { "192", range( 191.0, 193.0 ) },
    { "196", range( 195.0, 197.0 ) },
    { "200", range( 199.0, 201.0 ) },
    { "202", range( 201.0, 202.5 ) },
    { "205", range( 202.5, 205.5 ) },
    { "207", range( 205.5, 209.5 ) },
  } {
  SetBranchAddressChecked( "Il2mh", &nt_Il2mh ); 
  SetBranchAddressChecked( "Ebeam", &nt_Ebeam );
  SetBranchAddressChecked( "Pspr", &nt_Pspr );
  SetBranchAddressChecked( "Pspri", &nt_Pspri );
  SetBranchAddressChecked( "Lwqqln", &nt_Lwqqln );
  SetBranchAddressChecked( "Lwqqqq", &nt_Lwqqqq );
}

bool LEP2NtupleReader::Preselection( const std::string& ecms ) {
  bool result= true;
  Int_t icjst= nt_Icjst;
  Int_t iebst= nt_Iebst;
  Int_t il2mh= nt_Il2mh;
  if( icjst != 3 or iebst != 3 or il2mh != 1 ) result= false;
  range ecmsrange= ecmsranges.at( ecms );
  float ntecms= 2.0*nt_Ebeam;
  if( ( ntecms < ecmsrange.first ) or 
      ( ntecms > ecmsrange.second ) ) result= false;
  return result;
}

bool LEP2NtupleReader::Selection( const std::string& ecms ) {
  bool result= true;
  if( not Preselection( ecms ) ) result= false;
  if( nt_Ntkd02 < 7 ) result= false;
  if( abscostt() > 0.9 ) result= false;
  return result;
}

const std::map<std::string,bool> 
LEP2NtupleReader::getSelections( const std::string& ecms ) {
  std::map<std::string,bool> selections;
  selections["stand"]= false;
  selections["sprold"]= false;
  selections["costt07"]= false;
  if( not ( ecms == "130" or ecms == "136" ) ) {
    selections["wqqlnhi"]= false;
    selections["wqqlnlo"]= false;
    selections["wqqqqhi"]= false;
    selections["wqqqqlo"]= false;
  }
  if( Selection( ecms ) ) {
    double ecm= 2.0*nt_Ebeam;
    double sprime= sqrt( fabs( nt_Pspr[3]*nt_Pspr[3]-
			       nt_Pspr[0]*nt_Pspr[0]-
			       nt_Pspr[1]*nt_Pspr[1]-
			       nt_Pspr[2]*nt_Pspr[2] ) );
    double sprime_old= sqrt( fabs( nt_Pspri[3]*nt_Pspri[3]-
				   nt_Pspri[0]*nt_Pspri[0]-
				   nt_Pspri[1]*nt_Pspri[1]-
				   nt_Pspri[2]*nt_Pspri[2] ) );
    double sprcut= 10.0;
    bool lsprime= ecm - sprime < sprcut;
    bool lsprold= ecm - sprime_old < sprcut;
    if( ecms == "130" or ecms == "136" ) {
      if( lsprime ) selections["stand"]= true;
      if( lsprold ) selections["sprold"]= true;
      if( lsprime and abscostt() < 0.7 ) selections["costt07"]= true;
    }
    else {
      const double wqqln= 0.5;
      const double wqqlnHi= 0.75;
      const double wqqlnLo= 0.25;
      const double wqqqq= 0.25;
      const double wqqqqHi= 0.4;
      const double wqqqqLo= 0.1;
      bool lwqqln= nt_Lwqqln < wqqln;
      bool lwqqlnHi= nt_Lwqqln < wqqlnHi;
      bool lwqqlnLo= nt_Lwqqln < wqqlnLo;
      bool lwqqqq= nt_Lwqqqq < wqqqq;
      bool lwqqqqHi= nt_Lwqqqq < wqqqqHi;
      bool lwqqqqLo= nt_Lwqqqq < wqqqqLo;
      if( lsprime ) {
	if( lwqqln and lwqqqq ) selections["stand"]= true;
	if( lwqqln and lwqqqq and abscostt() < 0.7 ) selections["costt07"]= true;
	if( lwqqlnHi and lwqqqq ) selections["wqqlnhi"]= true;
	if( lwqqlnLo and lwqqqq ) selections["wqqlnlo"]= true;
	if( lwqqln and lwqqqqHi ) selections["wqqqqhi"]= true;
	if( lwqqln and lwqqqqLo ) selections["wqqqqlo"]= true;
      }
      if( lsprold and lwqqln and lwqqqq ) selections["sprold"]= true;
    }
  }
  return selections;
}
