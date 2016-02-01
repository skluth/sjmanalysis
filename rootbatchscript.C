{
gROOT->LoadMacro( "libNtupleReader.so");
gROOT->LoadMacro( "libNtupleReaderDict.so");
gROOT->ProcessLine(".include /home/skluth/qcd/fastjet/fastjet-3.1.3/install/include");
gROOT->LoadMacro( "LEP1Analysis.C+" );
LEP1Analysis(1000);
// LEP1Analysis();
}
