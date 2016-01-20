
#include "NtupleReader.hh"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
using std::cout;
using std::endl;
#include <sstream>
#include <string>
using std::string;
#include <stdexcept>
using std::logic_error;

// extern "C" {
//   void pxlth4_( Int_t*, Int_t*, Float_t*, Float_t*, Float_t*, Int_t* );
// };

NtupleReader::NtupleReader() : nt_file(0), nt_tree(0), nt_isMC(false), 
			       nt_vtlvcache(false) {}

NtupleReader::NtupleReader( const char* filename, const char* ntid ) :
  nt_file(0), nt_tree(0), nt_isMC(false), nt_vtlvcache(false) {
  OpenFileAndLoadNtuple( filename, ntid );
  return;
}

NtupleReader::~NtupleReader() {
  CloseFile();
}

void NtupleReader::OpenFileAndLoadNtuple( const char* filename, 
					  const char* ntid ) {
  cout << "NtupleReader::OpenFileAndLoadNtuple: opening file: " 
       << filename << endl;
  nt_file= new TFile( filename );
  if( not nt_file->IsOpen() ) {
    string txt= "NtupleReader::OpenFileAndLoadNtuple: file not open: ";
    txt+= filename;
    throw std::logic_error( txt );
  }
  nt_tree= (TTree*) nt_file->Get( ntid );
  if( nt_tree == 0 ) {
    string txt= "NtupleReader::OpenFileAndLoadNtuple: tree not found: ";
    txt+= ntid;
    throw std::logic_error( txt );
  }
  cout << "NtupleReader::OpenFileAndLoadNtuple: " 
       << GetNumberEntries() << " events on file" << endl;
  string sfilename( filename );
  if( sfilename.find( "mc" ) != string::npos ) nt_isMC= true;
  Init();
  nt_nevents= 0;
  return;
}

void NtupleReader::CloseFile() {
  if( nt_file ) {
    cout << "NtupleReader::CloseFile: closing file: " 
	 << nt_file->GetName() << ", " << nt_nevents << " events" << endl;
    nt_file->Close();
    delete nt_file;
    nt_file= 0;
  }
  else {
    throw std::logic_error( "NtupleReader::CloseFile: NULL file pointer" );
  }
  return;
}

Int_t NtupleReader::GetNumberEntries() {
  if( nt_tree ) return nt_tree->GetEntries();
  else return -1;
}

bool NtupleReader::GetEvent( Int_t ievnt ) {
  bool result= false;
  if( nt_tree and nt_tree->GetEvent( ievnt ) > 0 ) {
    result= true;
    nt_vtlvcache= false;
    nt_nevents++;
  }
  return result;
}

bool NtupleReader::LEP1Preselection() {
  bool result= true;
  Int_t icjst= nt_Icjst;
  Int_t iebst= nt_Iebst;
  Int_t itkmh= nt_Itkmh;
  if( icjst != 3 or iebst != 3 or itkmh != 1 ) result= false;
  return result;
}

bool NtupleReader::LEP1Selection() {
  bool result= true;
  if( !LEP1Preselection() ) result= false;
  if( nt_Ntkd02 < 5 ) result= false;
  if( costt() > 0.9 ) result= false;
  return result;
}

std::map<std::string,bool> NtupleReader::LEP1Selections() {
  std::map<std::string,bool> selections;
  selections["standard"]= false;
  selections["costt07"]= false;
  selections["nch7"]= false;
  bool preselection= LEP1Selection();
  if( preselection and nt_Ntkd02 >= 5 and costt() < 0.9 ) selections["stand"]= true;
  if( preselection and nt_Ntkd02 >= 7 and costt() < 0.9 ) selections["nch7"]= true;
  if( preselection and nt_Ntkd02 >= 5 and costt() < 0.7 ) selections["costt07"]= true;
  return selections;
}

bool NtupleReader::MCNonRad() {
  bool result= false;
  Int_t inonr= nt_Inonr;
  if( nt_isMC and inonr == 1 ) result= true;
  return result;
}

bool NtupleReader::inRange( Int_t njet, Int_t max ) {
  return njet > 0 and njet <= max;
}

Double_t NtupleReader::getRecoYmergeValue( const TString& reco, Int_t njet,
					   Int_t maxmt, Float_t* Ymt, 
					   Int_t maxtc, Float_t* Ytc, 
					   Int_t maxt, Float_t* Yt, 
					   Int_t maxc, Float_t* Yc, 
					   Int_t maxh, Float_t* Yh, 
					   Int_t maxp, Float_t* Yp ) {
  Double_t result= -1.0;
  if( reco == "mt" and inRange( njet, maxmt ) ) result= Ymt[njet-1];
  else if( reco == "tc" and inRange( njet, maxtc ) ) result= Ytc[njet-1];
  else if( reco == "tracks" and inRange( njet, maxt ) ) result= Yt[njet-1];
  else if( reco == "cluster" and inRange( njet, maxc ) ) result= Yc[njet-1];
  else if( reco == "hadron" and inRange( njet, maxh ) ) result= Yh[njet-1];
  else if( reco == "parton" and inRange( njet, maxp ) ) result= Yp[njet-1];
  else {
    std::cout << "NtupleReader::getRecoYmergeValue: no value found " << reco << "  " 
	      << njet << " " << maxmt << std::endl;
  }
  return result;
}

Double_t NtupleReader::getRecoValue( const TString& reco, 
				     Float_t mt,
				     Float_t tc,
				     Float_t tracks,
				     Float_t cluster,
				     Float_t hadron,
				     Float_t parton ) {
  Double_t value= -1.0;
  if( reco == "mt" ) value= mt;
  else if( reco == "tc" ) value= tc;
  else if( reco == "tracks" ) value= tracks;
  else if( reco == "clusters" ) value= cluster;
  else if( reco == "hadron" ) value= hadron;
  else if( reco == "parton" ) value= parton;
  else std::cout << "NtupleReader::getRecoValue: reco method " << reco << " not recognised" 
		 << std::endl;
  return value;
}

Double_t NtupleReader::getYmergeD( const TString& reco, Int_t njet ) {
  return getRecoYmergeValue( reco, njet, 
			     nt_Nxjdmt, nt_Yddmt, nt_Nxjdtc, nt_Yddtc,
			     nt_Nxjdt, nt_Yddt, nt_Nxjdc, nt_Yddc,
			     nt_Nxjdh, nt_Ydh, nt_Nxjdp, nt_Ydp );
}

Double_t NtupleReader::getYmergeE( const TString& reco, Int_t njet ) {
  return getRecoYmergeValue( reco, njet, 
			     nt_Nxjemt, nt_Yedmt, nt_Nxjetc, nt_Yedtc,
			     nt_Nxjet, nt_Yedt, nt_Nxjec, nt_Yedc,
			     nt_Nxjeh, nt_Yeh, nt_Nxjep, nt_Yep );
}

Double_t NtupleReader::getThrust( const TString& reco ) {
  return getRecoValue( reco, nt_Tdmt, nt_Tdtc, nt_Tdt, nt_Tdc, nt_Th, nt_Tp );
}

void NtupleReader::SetBranchAddressChecked( const char* branchname, void* address ) {
  if( nt_tree->GetBranch( branchname ) != 0 ) {
    nt_tree->SetBranchAddress( branchname, address );
  }
  else {
    std::cout << "Can't set branch address " << branchname << std::endl;
  }
  return;
}

void NtupleReader::Init() {

  // LEP1 preselection:
  nt_tree->SetBranchAddress( "Icjst", &nt_Icjst );
  nt_tree->SetBranchAddress( "Iebst", &nt_Iebst );
  nt_tree->SetBranchAddress( "Itkmh", &nt_Itkmh );

  // LEP1 selection:
  nt_tree->SetBranchAddress( "Ntkd02", &nt_Ntkd02 );
  nt_tree->SetBranchAddress( "Tvectc", &nt_Tvectc );

  // Calculated shapes and jets:
  nt_tree->SetBranchAddress( "Tdmt", &nt_Tdmt );
  nt_tree->SetBranchAddress( "Tdt", &nt_Tdt );
  nt_tree->SetBranchAddress( "Tdc", &nt_Tdc );
  nt_tree->SetBranchAddress( "Tdtc", &nt_Tdtc );
  nt_tree->SetBranchAddress( "Nxjdmt", &nt_Nxjdmt );
  nt_tree->SetBranchAddress( "Yddmt", &nt_Yddmt );
  nt_tree->SetBranchAddress( "Nxjdt", &nt_Nxjdt );
  nt_tree->SetBranchAddress( "Yddt", &nt_Yddt );
  nt_tree->SetBranchAddress( "Nxjdc", &nt_Nxjdc );
  nt_tree->SetBranchAddress( "Yddc", &nt_Yddc );
  nt_tree->SetBranchAddress( "Nxjdtc", &nt_Nxjdtc );
  nt_tree->SetBranchAddress( "Yddtc", &nt_Yddtc );

  nt_tree->SetBranchAddress( "Nxjemt", &nt_Nxjemt );
  nt_tree->SetBranchAddress( "Yedmt", &nt_Yedmt );
  nt_tree->SetBranchAddress( "Nxjet", &nt_Nxjet );
  nt_tree->SetBranchAddress( "Yedt", &nt_Yedt );
  nt_tree->SetBranchAddress( "Nxjec", &nt_Nxjec );
  nt_tree->SetBranchAddress( "Yedc", &nt_Yedc );
  nt_tree->SetBranchAddress( "Nxjetc", &nt_Nxjetc );
  nt_tree->SetBranchAddress( "Yedtc", &nt_Yedtc );

  // Tracks and clusters:
  nt_tree->SetBranchAddress( "Ntrk", &nt_Ntrk );
  nt_tree->SetBranchAddress( "Id02", &nt_Id02 );
  nt_tree->SetBranchAddress( "Ptrk", &nt_Ptrk );
  nt_tree->SetBranchAddress( "Nclus", &nt_Nclus );
  nt_tree->SetBranchAddress( "Pclus", &nt_Pclus );
  nt_tree->SetBranchAddress( "Nmttrk", &nt_Nmttrk );
  nt_tree->SetBranchAddress( "Imttrk", &nt_Imttrk );
  nt_tree->SetBranchAddress( "Mtscft", &nt_Mtscft );
  nt_tree->SetBranchAddress( "Nmtcls", &nt_Nmtcls );
  nt_tree->SetBranchAddress( "Nmtkil", &nt_Nmtkil );
  nt_tree->SetBranchAddress( "Imtkil", &nt_Imtkil );
  nt_tree->SetBranchAddress( "Imtcls", &nt_Imtcls );
  nt_tree->SetBranchAddress( "Mtscfc", &nt_Mtscfc );

  // MC quantities:

  // Partons and hadrons:
  if( nt_isMC ) {
    SetBranchAddressChecked( "Inonr", &nt_Inonr );
    SetBranchAddressChecked( "Ntrkp", &nt_Ntrkp );
    SetBranchAddressChecked( "Ptrkp", &nt_Ptrkp );
    SetBranchAddressChecked( "Ntrkh", &nt_Ntrkh );
    SetBranchAddressChecked( "Ptrkh", &nt_Ptrkh );

    nt_tree->SetBranchAddress( "Th", &nt_Th );
    nt_tree->SetBranchAddress( "Tp", &nt_Tp );
    nt_tree->SetBranchAddress( "Nxjdh", &nt_Nxjdh );
    nt_tree->SetBranchAddress( "Ydh", &nt_Ydh );
    nt_tree->SetBranchAddress( "Nxjdp", &nt_Nxjdp );
    nt_tree->SetBranchAddress( "Ydp", &nt_Ydp );
    nt_tree->SetBranchAddress( "Nxjeh", &nt_Nxjeh );
    nt_tree->SetBranchAddress( "Yeh", &nt_Yeh );
    nt_tree->SetBranchAddress( "Nxjep", &nt_Nxjep );
    nt_tree->SetBranchAddress( "Yep", &nt_Yep );

  }

  return;

}  

const std::vector<TLorentzVector>& NtupleReader::GetLorentzVectors( const std::string & opt ) {
  static std::vector<TLorentzVector> vtlv;
  static std::string lastopt= "none";
  if( nt_vtlvcache and lastopt == opt ) return vtlv;
  static Float_t ptrack[nt_maxtrk][4];
  Int_t ntrack;
  if( opt == "parton" ) {
    GetP( ptrack, nt_maxtrk, ntrack );
  }
  else if( opt == "hadron" ) {
    GetH( ptrack, nt_maxtrk, ntrack );
  }
  else if( opt == "tracks" ) {
    GetTrk( ptrack, nt_maxtrk, ntrack );
  }
  else if( opt == "clusters" ) {
    GetCls( ptrack, nt_maxtrk, ntrack );
  }
  else if( opt == "tc" ) {
    GetTC( ptrack, nt_maxtrk, ntrack );
  }
  else if( opt == "mt" ) {
    GetMt( ptrack, nt_maxtrk, ntrack );
  }
  else {
    std::cout << "NtupleReader::GetLorentzVectors: option " << opt << " not recognised" 
	      << std::endl;
    exit( 1 );
  }
  vtlv.clear();
  vtlv.reserve( ntrack );
  for( Int_t itrk= 0; itrk < ntrack; itrk++ ) {
    TLorentzVector tlv( ptrack[itrk][0], ptrack[itrk][1], 
			ptrack[itrk][2], ptrack[itrk][3] );
    vtlv.push_back( tlv );
  }
  nt_vtlvcache= true;
  lastopt= opt;
  return vtlv;
}

void NtupleReader::GetP( Float_t ptrack[][4], Int_t maxtrack, Int_t & ntrack ) {
  for( Int_t itrk= 0; itrk < nt_Ntrkp; itrk++ ) {
    for( Int_t j= 0; j < 4; j++ ) {
      ptrack[itrk][j]= nt_Ptrkp[itrk][j];
    }
  }
  ntrack= nt_Ntrkp;
  return;
}
void NtupleReader::GetH( Float_t ptrack[][4], Int_t maxtrack, Int_t & ntrack ) {
  for( Int_t itrk= 0; itrk < nt_Ntrkh; itrk++ ) {
    for( Int_t j= 0; j < 4; j++ ) {
      ptrack[itrk][j]= nt_Ptrkh[itrk][j];
    }
  }
  ntrack= nt_Ntrkh;
  return;
}

void NtupleReader::GetTC( Float_t ptrack[][4], Int_t maxtrack, Int_t & ntrack ) {
  GetTrk( ptrack, maxtrack, ntrack );
  GetCls( ptrack, maxtrack, ntrack, ntrack );
  return;
}

void NtupleReader::GetTrk( Float_t ptrack[][4], Int_t maxtrack, Int_t & ntrack ) {
  Float_t mpi2= pow( 0.140, 2 );
  Int_t ifill= 0;
  for( Int_t itrk= 0; itrk < nt_Ntrk; itrk++ ) {
    if( nt_Id02[itrk] == 0 ) continue;
    if( ifill == maxtrack ) {
      std::cout << "NtupleReader::getTrk: array too small " << ifill << std::endl;
      break;
    }
    Float_t sum= 0.0;
    for( Int_t j= 0; j < 3; j++ ) {
      ptrack[ifill][j]= nt_Ptrk[itrk][j];
      sum+= pow( ptrack[ifill][j], 2 );
    }
    ptrack[ifill][3]= sqrt( sum+mpi2 );
    ifill++;
  }
  ntrack= ifill;
  return;
}
void NtupleReader::GetCls( Float_t ptrack[][4], Int_t maxtrack, Int_t & ntrack,
			   Int_t ioff ) {
  if( nt_Nclus > maxtrack ) {
    std::cout << "Ntuple::getcls: array too small " << maxtrack << std::endl;
  }
  Int_t iclus;
  for( iclus= 0; iclus < TMath::Min( maxtrack, nt_Nclus ); iclus++ ) {
    Float_t sum= 0.0;
    for( Int_t j= 0; j < 3; j++ ) {
      ptrack[ioff+iclus][j]= nt_Pclus[iclus][j];
      sum+= pow( ptrack[ioff+iclus][j], 2 );
    }
    ptrack[ioff+iclus][3]= sqrt( sum );
  }
  ntrack= ioff+iclus;
  return;
}
void NtupleReader::GetMt( Float_t ptrack[][4], Int_t maxtrack, Int_t & ntrack ) {

  // Tracks first:
  Float_t mpi2= pow( 0.140, 2 );
  Int_t ifill= 0;
  for( Int_t itrk= 0; itrk < nt_Ntrk; itrk++ ) {
    if( ifill == maxtrack ) {
      std::cout << "NtupleReader::getMt: array too small " << ifill << std::endl;
      break;
    }
    // Check if track is selected:
    if( nt_Id02[itrk] == 0 ) continue;
    // Check if track is scaled:
    Float_t scf= 1.0;
    for( Int_t jmttrk= 0; jmttrk < nt_Nmttrk; jmttrk++ ) {
      if( nt_Imttrk[jmttrk]-1 == itrk ) {
	scf= nt_Mtscft[jmttrk];
	break;
      }
    }
    // Copy track components:
    Float_t sum= 0.0; 
    for( Int_t j= 0; j < 3; j++ ) {
      ptrack[ifill][j]= nt_Ptrk[itrk][j]*scf;
      sum+= pow( ptrack[ifill][j], 2 );
    }
    ptrack[ifill][3]= sqrt( sum + mpi2 );
    ifill++;
  }    
 
  // Clusters are either killed, scaled or copied:
  for( Int_t iclus= 0; iclus < TMath::Min( maxtrack, nt_Nclus ); iclus++ ) {
    if( ifill == maxtrack ) {
      std::cout << "Ntuple::getmt: array too small " << ifill << std::endl;
      break;
    }
    // Check if cluster is killed:
    bool killed= false;
    for( Int_t jmtkil= 0; jmtkil < nt_Nmtkil; jmtkil++ ) {
      if( nt_Imtkil[jmtkil]-1 == iclus ) {
	killed= true;
	break;
      }
    }
    if( killed ) continue;
    // Check if cluster is scaled:
    Float_t scf= 1.0;
    for( Int_t jmtcls= 0; jmtcls < nt_Nmtcls; jmtcls++ ) {
      if( nt_Imtcls[jmtcls]-1 == iclus ) {
	scf= nt_Mtscfc[jmtcls];
	break;
      }
    }
    // Copy cluster components:
    Float_t sum= 0.0;
    for( Int_t j= 0; j < 3; j++ ) {
      ptrack[ifill][j]= nt_Pclus[iclus][j]*scf;
      sum+= pow( ptrack[ifill][j], 2 );
    }
    ptrack[ifill][3]= sqrt( sum );
    ifill++;

  }

  // The End:
  ntrack= ifill;
  return;

}

