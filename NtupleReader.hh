#ifndef NTUPLEREADER_HH
#define NTUPLEREADER_HH

#include "Rtypes.h"
#include "TLorentzVector.h"
#include "TString.h"
#include <vector>
#include <string>
#include <map>

class TFile;
class TTree;

class NtupleReader {

public:

  NtupleReader();
  NtupleReader( const char* filename, const char* ntid="h10");
  ~NtupleReader();

  bool OpenFileAndLoadNtuple( const char* filename, const char* ntid="h10" );
  void CloseFile();

  Int_t GetNumberEntries();
  bool GetEvent( Int_t ievnt );
 
  const std::vector<TLorentzVector>& GetLorentzVectors( const std::string & opt );
  void GetP( Float_t ptrack[][4], Int_t maxtrack, Int_t & ntrack );
  void GetH( Float_t ptrack[][4], Int_t maxtrack, Int_t & ntrack );
  void GetTrk( Float_t ptrack[][4], Int_t maxtrack, Int_t & ntrack );
  void GetCls( Float_t ptrack[][4], Int_t maxtrack, Int_t & ntrack, Int_t ioff=0 );
  void GetTC( Float_t ptrack[][4], Int_t maxtrack, Int_t & ntrack );
  void GetMt( Float_t ptrack[][4], Int_t maxtrack, Int_t & ntrack );

  bool LEP1Preselection();
  bool LEP1Selection();
  std::map<std::string,bool> LEP1Selections();
  bool MCNonRad();
  bool isMC() { return nt_isMC; }

  Float_t costt() { return nt_Tvectc[2]; }
  Double_t getYmergeD( const TString& reco, Int_t njet );
  Double_t getYmergeE( const TString& reco, Int_t njet );
  Double_t getThrust( const TString& reco );

private:

  void Init();
  void SetBranchAddressChecked( const char*, void* );
  bool inRange( Int_t, Int_t );
  Double_t getRecoValue( const TString&, 
			 Float_t, Float_t, Float_t, Float_t, Float_t, Float_t );
  Double_t getRecoYmergeValue( const TString& reco, Int_t njet,
			       Int_t maxmt, Float_t* Ymt, 
			       Int_t maxtc, Float_t* Ytc, 
			       Int_t maxt, Float_t* Yt, 
			       Int_t maxc, Float_t* Yc, 
			       Int_t maxh, Float_t* Yh, 
			       Int_t maxp, Float_t* Yp );

  TFile* nt_file;
  TTree* nt_tree;
  bool nt_isMC;
  bool nt_vtlvcache;

  static const Int_t nt_maxtrk= 501;
  static const Int_t nt_maxp= 50;
  static const Int_t nt_maxh= 2004;
  UShort_t nt_Irun;
  Int_t    nt_Ievnt;
  UChar_t  nt_Itkmh;
  UChar_t  nt_Igpmh;
  UChar_t  nt_Isist;
  UChar_t  nt_Icvst;
  UChar_t  nt_Icjst;
  UChar_t  nt_Iczst;
  UChar_t  nt_Iebst;
  UChar_t  nt_Ieest;
  UChar_t  nt_Istg1;
  UShort_t nt_Ntkd02;
  Float_t  nt_Ebeam;
  Float_t  nt_Pgce[4];
  Float_t  nt_Tvectc[3];
  UChar_t  nt_Il2mh;
  UChar_t  nt_Iwqqqq;
  Char_t   nt_Iwqqln;
  Int_t    nt_Idqqln;
  UChar_t  nt_Iebal;
  UChar_t  nt_Icat;
  UChar_t  nt_Iminrmp;
  Float_t  nt_Pspri[4];
  Float_t  nt_Pspr[4];
  Float_t  nt_Ephot;
  Float_t  nt_Efit;
  Float_t  nt_Pfit;
  Float_t  nt_Efd;
  Float_t  nt_W310;
  Float_t  nt_W420;
  Float_t  nt_Lwqqqq;
  Float_t  nt_Lwqqln;
  Float_t  nt_Lwexef0;
  Float_t  nt_Lwexef1;
  Float_t  nt_Tdtc;
  Float_t  nt_Tmadtc;
  Float_t  nt_Tmidtc;
  Float_t  nt_Mhdtc;
  Float_t  nt_Mldtc;
  Float_t  nt_Btdtc;
  Float_t  nt_Bwdtc;
  Float_t  nt_Cpdtc;
  Float_t  nt_Dpdtc;
  Float_t  nt_Sdtc;
  Float_t  nt_Adtc;
  Float_t  nt_Acpdtc;
  Float_t  nt_Tdt;
  Float_t  nt_Tmadt;
  Float_t  nt_Tmidt;
  Float_t  nt_Mhdt;
  Float_t  nt_Mldt;
  Float_t  nt_Btdt;
  Float_t  nt_Bwdt;
  Float_t  nt_Cpdt;
  Float_t  nt_Dpdt;
  Float_t  nt_Sdt;
  Float_t  nt_Adt;
  Float_t  nt_Acpdt;
  Float_t  nt_Tdc;
  Float_t  nt_Tmadc;
  Float_t  nt_Tmidc;
  Float_t  nt_Mhdc;
  Float_t  nt_Mldc;
  Float_t  nt_Btdc;
  Float_t  nt_Bwdc;
  Float_t  nt_Cpdc;
  Float_t  nt_Dpdc;
  Float_t  nt_Sdc;
  Float_t  nt_Adc;
  Float_t  nt_Acpdc;
  Float_t  nt_Tdmt;
  Float_t  nt_Tmadmt;
  Float_t  nt_Tmidmt;
  Float_t  nt_Mhdmt;
  Float_t  nt_Mldmt;
  Float_t  nt_Btdmt;
  Float_t  nt_Bwdmt;
  Float_t  nt_Cpdmt;
  Float_t  nt_Dpdmt;
  Float_t  nt_Sdmt;
  Float_t  nt_Admt;
  Float_t  nt_Acpdmt;
  Int_t    nt_Nxjdtc;
  Int_t    nt_Nxjdt;
  Int_t    nt_Nxjdc;
  Int_t    nt_Nxjdmt;
  Int_t    nt_Nxjetc;
  Int_t    nt_Nxjet;
  Int_t    nt_Nxjec;
  Int_t    nt_Nxjemt;
  Int_t    nt_Nxjctc;
  Int_t    nt_Nxjct;
  Int_t    nt_Nxjcc;
  Int_t    nt_Nxjcmt;
  Float_t  nt_Yddtc[31];   //[Nxjdtc]
  Float_t  nt_Yedtc[31];   //[Nxjetc]
  Float_t  nt_Ycdtc[31];   //[Nxjctc]
  Char_t   nt_Njcedtc[7];
  Char_t   nt_Njcrdtc[7];
  Float_t  nt_Yddt[31];   //[Nxjdt]
  Float_t  nt_Yedt[31];   //[Nxjet]
  Float_t  nt_Ycdt[31];   //[Nxjct]
  Char_t   nt_Njcedt[7];
  Char_t   nt_Njcrdt[7];
  Float_t  nt_Yddc[31];   //[Nxjdc]
  Float_t  nt_Yedc[31];   //[Nxjec]
  Float_t  nt_Ycdc[31];   //[Nxjcc]
  Char_t   nt_Njcedc[7];
  Char_t   nt_Njcrdc[7];
  Float_t  nt_Yddmt[31];   //[Nxjdmt]
  Float_t  nt_Yedmt[31];   //[Nxjemt]
  Float_t  nt_Ycdmt[31];   //[Nxjcmt]
  Char_t   nt_Njcedmt[7];
  Char_t   nt_Njcrdmt[7];
  Int_t    nt_Ntrk;
  UChar_t  nt_Id02[nt_maxtrk];   //[Ntrk]
  Float_t  nt_Dedx[nt_maxtrk];   //[Ntrk]
  Float_t  nt_Dded[nt_maxtrk];   //[Ntrk]
  UChar_t  nt_Nhde[nt_maxtrk];   //[Ntrk]
  Float_t  nt_Dp[nt_maxtrk];   //[Ntrk]
  Float_t  nt_Ptrk[nt_maxtrk][3];   //[Ntrk]
  Char_t   nt_Ichg[nt_maxtrk];   //[Ntrk]
  UChar_t  nt_Nhcj[nt_maxtrk];   //[Ntrk]
  Float_t  nt_Z0[nt_maxtrk];   //[Ntrk]
  Float_t  nt_D0[nt_maxtrk];   //[Ntrk]
  Int_t    nt_Nmttrk;
  UShort_t nt_Imttrk[nt_maxtrk];   //[Nmttrk]
  Float_t  nt_Mtscft[nt_maxtrk];   //[Nmttrk]
  Int_t    nt_Nclus;
  Int_t    nt_Nmtcls;
  UShort_t nt_Imtcls[nt_maxtrk];   //[Nmtcls]
  Int_t    nt_Nmtkil;
  UShort_t nt_Imtkil[nt_maxtrk];   //[Nmtkil]
  Float_t  nt_Pclus[1503][3];   //[Nclus]
  Float_t  nt_Mtscfc[nt_maxtrk];   //[Nmtcls]
  UChar_t  nt_Ioselbt;
  UChar_t  nt_Levslbt;
  Int_t    nt_Nvtxbt;
  UChar_t  nt_Ivmulbt[30];   //[Nvtxbt]
  UChar_t  nt_Nvsigbt[30];   //[Nvtxbt]
  Float_t  nt_Thrvecbt[3];
  Float_t  nt_Prvtxbt[3];
  Float_t  nt_Vnnbt[30];   //[Nvtxbt]
  Float_t  nt_Vchi2bt[30];   //[Nvtxbt]
  Float_t  nt_Vtxbt[30][3];   //[Nvtxbt]
  Float_t  nt_Pvtxbt[30][5];   //[Nvtxbt]
  Float_t  nt_Vdlen3bt[30];   //[Nvtxbt]
  Float_t  nt_Vderr3bt[30];   //[Nvtxbt]
  UChar_t  nt_Ievtyp;
  UChar_t  nt_Inonr;
  Float_t  nt_Pisr[4];
  Int_t    nt_Nprimf;
  Char_t   nt_Iferid[4];   //[Nprimf]
  Float_t  nt_Primf[4][4];   //[Nprimf]
  Float_t  nt_Tp;
  Float_t  nt_Tmap;
  Float_t  nt_Tmip;
  Float_t  nt_Mhp;
  Float_t  nt_Mlp;
  Float_t  nt_Btp;
  Float_t  nt_Bwp;
  Float_t  nt_Cpp;
  Float_t  nt_Dpp;
  Float_t  nt_Sp;
  Float_t  nt_Ap;
  Float_t  nt_Acpp;
  Float_t  nt_Th;
  Float_t  nt_Tmah;
  Float_t  nt_Tmih;
  Float_t  nt_Mhh;
  Float_t  nt_Mlh;
  Float_t  nt_Bth;
  Float_t  nt_Bwh;
  Float_t  nt_Cph;
  Float_t  nt_Dph;
  Float_t  nt_Sh;
  Float_t  nt_Ah;
  Float_t  nt_Acph;
  Int_t    nt_Nxjdp;
  Int_t    nt_Nxjdh;
  Int_t    nt_Nxjep;
  Int_t    nt_Nxjeh;
  Int_t    nt_Nxjcp;
  Int_t    nt_Nxjch;
  Float_t  nt_Ydp[31];   //[Nxjdp]
  Float_t  nt_Yep[31];   //[Nxjep]
  Float_t  nt_Ycp[31];   //[Nxjcp]
  Char_t   nt_Njcep[7];
  Char_t   nt_Njcrp[7];
  Float_t  nt_Ydh[31];   //[Nxjdh]
  Float_t  nt_Yeh[31];   //[Nxjeh]
  Float_t  nt_Ych[31];   //[Nxjch]
  Char_t   nt_Njceh[7];
  Char_t   nt_Njcrh[7];
  Int_t    nt_Ntrkp;
  Char_t   nt_Ilucp[nt_maxp];   //[Ntrkp]
  Float_t  nt_Ptrkp[nt_maxp][4];   //[Ntrkp]
  Int_t    nt_Ntrkh;
  Int_t    nt_Ntrk2;
  Float_t  nt_Ptrkh[nt_maxh][4];   //[Ntrkh]
  Int_t    nt_Iluch[nt_maxh];   //[Ntrkh]
  Int_t    nt_Iluc[nt_maxtrk];   //[Ntrk2]
  UChar_t  nt_Istrt[nt_maxtrk];   //[Ntrk2]
  Char_t   nt_Ichgh[nt_maxh];   //[Ntrkh]
  
};

#endif
