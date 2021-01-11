
#include "LEPNtupleReader.hh"
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
#include <numeric>

// Link PXLIB Fortran
// extern "C" {
//   void pxlth4_( Int_t*, Int_t*, Float_t*, Float_t*, Float_t*, Int_t* );
// };

LEPNtupleReader::LEPNtupleReader() : nt_file(0), nt_tree(0), nt_isMC(false), lprint(false) {}

LEPNtupleReader::LEPNtupleReader( const char* filename, const char* ntid, const bool lpr ) :
  nt_file(0), nt_tree(0), nt_isMC(false), nt_vtlvcache{false},
  nt_nevents(0), lprint(lpr),
  vtlvCache{ { "parton", std::vector<TLorentzVector>() },
    { "hadron", std::vector<TLorentzVector>() },
    { "tracks", std::vector<TLorentzVector>() },
    { "cluster", std::vector<TLorentzVector>() },
    { "tc", std::vector<TLorentzVector>() },
    { "mt", std::vector<TLorentzVector>() } },
  cacheIsValid{ { "parton", false },
    { "hadron", false },
    { "tracks", false },
    { "cluster", false },
    { "tc", false },
    { "mt", false } }
{
  OpenFileAndLoadNtuple( filename, ntid );
  return;
}

LEPNtupleReader::~LEPNtupleReader() {
  try { CloseFile(); }
  catch( const std::runtime_error& e ) {}
}

void LEPNtupleReader::OpenFileAndLoadNtuple( const char* filename, 
					     const char* ntid ) {
  if( lprint ) {
    cout << "LEPNtupleReader::OpenFileAndLoadNtuple: opening file: " 
	 << filename << endl;
  }
  nt_file= TFile::Open( filename );
  if( not nt_file ) {
    string txt= "LEPNtupleReader::OpenFileAndLoadNtuple: file not found: ";
    txt+= filename;
    throw std::runtime_error( txt );
  }  
  if( not nt_file->IsOpen() ) {
    string txt= "LEPNtupleReader::OpenFileAndLoadNtuple: file not open: ";
    txt+= filename;
    throw std::runtime_error( txt );
  }
  nt_tree= (TTree*) nt_file->Get( ntid );
  if( nt_tree == 0 ) {
    string txt= "LEPNtupleReader::OpenFileAndLoadNtuple: tree not found: ";
    txt+= ntid;
    throw std::runtime_error( txt );
  }
  if( lprint ) {
    cout << "LEPNtupleReader::OpenFileAndLoadNtuple: " 
	 << GetNumberEntries() << " events on file" << endl;
  }
  string sfilename( filename );
  if( sfilename.find( "mc" ) != string::npos ) nt_isMC= true;
  Init();
  nt_nevents= 0;
  return;
}

void LEPNtupleReader::CloseFile() {
  if( nt_file ) {
    if( lprint ) {
      cout << "LEPNtupleReader::CloseFile: closing file: " 
	   << nt_file->GetName() << ", " << nt_nevents << " events" << endl;
    }
    nt_file->Close();
    delete nt_file;
    nt_file= 0;
  }
  else {
    throw std::runtime_error( "LEPNtupleReader::CloseFile: NULL file pointer" );
  }
  return;
}

Int_t LEPNtupleReader::GetNumberEntries() {
  if( nt_tree ) return nt_tree->GetEntries();
  else return -1;
}

bool LEPNtupleReader::GetNextEvent( Int_t maxevt ) {
  bool maxreached= maxevt > 0 and maxevt == nt_nevents;
  return not maxreached and GetEvent( nt_nevents );
}

bool LEPNtupleReader::GetEvent( Int_t ievent ) {
  bool result= false;
  if( nt_tree and nt_tree->GetEvent( ievent ) > 0 ) {
    result= true;
    for( const auto & keyValue : cacheIsValid ) cacheIsValid[keyValue.first]= false;
    nt_nevents++;
  }
  return result;
}

bool LEPNtupleReader::MCNonRad() {
  bool result= false;
  Int_t inonr= nt_Inonr;
  if( nt_isMC and inonr == 1 ) result= true;
  return result;
}

bool LEPNtupleReader::inRange( Int_t njet, Int_t max ) {
  return njet > 0 and njet <= max;
}

Double_t LEPNtupleReader::getRecoYmergeValue( const std::string& reco, 
					      Int_t njet,
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
  // else {
  //   std::cout << "LEPNtupleReader::getRecoYmergeValue: no value found " << reco << " "
  // 	      << njet << " ";
  //   if( reco == "mt" ) std::cout << maxmt;
  //   else if( reco == "tc" ) std::cout << maxtc;
  //   else if( reco == "tracks" ) std::cout << maxt;
  //   else if( reco == "cluster" ) std::cout << maxc;
  //   else if( reco == "hadron" ) std::cout << maxh;
  //   else if( reco == "parton" ) std::cout << maxp;
  //   std::cout << std::endl;
  // }
  return result;
}

Double_t LEPNtupleReader::getRecoValue( const std::string& reco, 
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
  else std::cout << "LEPNtupleReader::getRecoValue: reco method " 
		 << reco << " not recognised" << std::endl;
  return value;
}

Double_t LEPNtupleReader::getYmerge( const std::string& algorithm, 
				     const std::string& reco, Int_t njet ) {
  if( algorithm == "jade" ) {
    return getRecoYmergeValue( reco, njet, 
			       nt_Nxjemt, nt_Yedmt, nt_Nxjetc, nt_Yedtc,
			       nt_Nxjet, nt_Yedt, nt_Nxjec, nt_Yedc,
			       nt_Nxjeh, nt_Yeh, nt_Nxjep, nt_Yep );
  }
  else if( algorithm == "durham" ) {
    return getRecoYmergeValue( reco, njet, 
			       nt_Nxjdmt, nt_Yddmt, nt_Nxjdtc, nt_Yddtc,
			       nt_Nxjdt, nt_Yddt, nt_Nxjdc, nt_Yddc,
			       nt_Nxjdh, nt_Ydh, nt_Nxjdp, nt_Ydp );    
  }
  else {
    throw std::runtime_error( "LEPNtupleReader::getYmerge: algorithm not known "+algorithm );
  }
}


Double_t LEPNtupleReader::getThrust( const std::string& reco ) {
  return getRecoValue( reco, nt_Tdmt, nt_Tdtc, nt_Tdt, nt_Tdc, nt_Th, nt_Tp );
}

void LEPNtupleReader::SetBranchAddressChecked( const char* branchname, void* address ) {
  if( nt_tree->GetBranch( branchname ) != 0 ) {
    nt_tree->SetBranchAddress( branchname, address );
  }
  else {
    std::cout << "Can't set branch address " << branchname << std::endl;
  }
  return;
}

void LEPNtupleReader::Init() {

  // Preselection:
  nt_tree->SetBranchAddress( "Icjst", &nt_Icjst );
  nt_tree->SetBranchAddress( "Iebst", &nt_Iebst );
  // ITKM and L2MH in subclasses
  //  nt_tree->SetBranchAddress( "Itkmh", &nt_Itkmh );

  // General selection:
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

const std::vector<TLorentzVector> LEPNtupleReader::GetLorentzVectors( const std::string & opt ) {
  if( cacheIsValid[opt] ) return vtlvCache[opt];
  vtlv= vtlvCache[opt];
  if( opt == "parton" ) {
    getPTlv();
  }
  else if( opt == "hadron" ) {
    getHTlv();
  }
  else if( opt == "tracks" ) {
    getTrkTlv();
  }
  else if( opt == "clusters" ) {
    getClsTlv();
  }
  else if( opt == "tc" ) {
    getTCTlv();
  }
  else if( opt == "mt" ) {
    getMtTlv();
  }
  else {
    throw std::runtime_error( "LEPNtupleReader::GetLorentzVectors: option "+opt+" not recognised" );
  }
  vtlvCache[opt]= vtlv;
  cacheIsValid[opt]= true;
  return vtlv;
}

void LEPNtupleReader::getPTlv() {
  vtlv.resize( nt_Ntrkp );
  for( Int_t itrk= 0; itrk < nt_Ntrkp; itrk++ ) {
    for( Int_t j= 0; j < 4; j++ ) vtlv[itrk][j]= nt_Ptrkp[itrk][j];
  }
}

void LEPNtupleReader::getHTlv() {
  vtlv.resize( nt_Ntrkh );
  for( Int_t itrk= 0; itrk < nt_Ntrkh; itrk++ ) {
    for( Int_t j= 0; j < 4; j++ ) vtlv[itrk][j]= nt_Ptrkh[itrk][j];
  }
}

void LEPNtupleReader::getTCTlv() {
  Int_t ntrack= getTrkTlv();
  getClsTlv( ntrack );
  return;
}

Int_t LEPNtupleReader::getTrkTlv() {
  vtlv.resize( nt_Ntrk );
  Float_t mpi2= pow( 0.140, 2 );
  Int_t ifill= 0;
  for( Int_t itrk= 0; itrk < nt_Ntrk; itrk++ ) {
    if( nt_Id02[itrk] == 0 ) continue;
    if( ifill == nt_maxtrk ) {
      std::cout << "LEPNtupleReader::getTrk: array too small " << ifill << std::endl;
      break;
    }
    Float_t sum= 0.0;
    for( Int_t j= 0; j < 3; j++ ) {
      vtlv[ifill][j]= nt_Ptrk[itrk][j];
      sum+= pow( vtlv[ifill][j], 2 );
    }
    vtlv[ifill][3]= sqrt( sum+mpi2 );
    ifill++;
  }
  vtlv.resize( ifill );
  return ifill;
}

void LEPNtupleReader::getClsTlv( Int_t ioff ) {
  vtlv.resize( ioff+nt_Nclus );
  for( Int_t iclus= 0; iclus < nt_Nclus; iclus++ ) {
    Float_t sum= 0.0;
    for( Int_t j= 0; j < 3; j++ ) {
      vtlv[ioff+iclus][j]= nt_Pclus[iclus][j];
      sum+= pow( vtlv[ioff+iclus][j], 2 );
    }
    vtlv[ioff+iclus][3]= sqrt( sum );
  }
  return;
}

void LEPNtupleReader::getMtTlv() {

  vtlv.resize( nt_Ntrk+nt_Nclus );

  // Tracks first:
  Float_t mpi2= pow( 0.140, 2 );
  Int_t ifill= 0;
  for( Int_t itrk= 0; itrk < nt_Ntrk; itrk++ ) {
    if( ifill == nt_maxtrk ) {
      std::cout << "LEPNtupleReader::getMtTlv: array too small " << ifill << std::endl;
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
      vtlv[ifill][j]= nt_Ptrk[itrk][j]*scf;
      sum+= pow( vtlv[ifill][j], 2 );
    }
    vtlv[ifill][3]= sqrt( sum + mpi2 );
    ifill++;
  }

  // Clusters are either killed, scaled or copied:
  for( Int_t iclus= 0; iclus < TMath::Min( nt_maxtrk, nt_Nclus ); iclus++ ) {
    if( ifill == nt_maxtrk ) {
      std::cout << "Ntuple::getMtTlv: array too small " << ifill << std::endl;
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
      vtlv[ifill][j]= nt_Pclus[iclus][j]*scf;
      sum+= pow( vtlv[ifill][j], 2 );
    }
    vtlv[ifill][3]= sqrt( sum );
    ifill++;

  }

  // The End:
  vtlv.resize( ifill );
  return;

}


Double_t LEPNtupleReader::Evis( const std::vector<TLorentzVector>& v ) const {
  Double_t evis= std::accumulate( v.begin(), v.end(), 0.0,
				  []( Double_t sum, const TLorentzVector& tlv ) {
				    return sum+= tlv.E();
				  }
				  );
  return evis;
}
