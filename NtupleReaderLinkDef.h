//#ifdef __MAKECINT__
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class vector<TLorentzVector>+;
#pragma link C++ class vector<string>+;
#pragma link C++ class map<string,Bool_t>+;
#pragma link C++ class pair<string,bool>+;
#pragma link C++ class map<TString,Int_t>+;
#pragma link C++ class vector<Observable*>+;
#pragma link C++ class vector<FilledObservable*>+;
#pragma link C++ class map<string,DataStructure*>+;
#pragma link C++ class map<string,MatrixDataStructure*>+;
#pragma link C++ class vector<Analysis>+;
#pragma link C++ class NtupleReader;
#pragma link C++ class Analysis;
#pragma link C++ class DataStructure;
#pragma link C++ class MatrixDataStructure;
#pragma link C++ class JetrateDataStructure;
#pragma link C++ class DifferentialDataStructure;
#pragma link C++ class Observable;
#pragma link C++ class FilledObservable;
#pragma link C++ class ObservableFactory;
#pragma link C++ class ObsFastJetDiff;
#pragma link C++ class ObsPartonShower;
#pragma link C++ class ObsDifferential;
#pragma link C++ class ObsJetrate;
#pragma link C++ class TFastJet;
#pragma link C++ class Unfolder;
#pragma link C++ class OutputWriter;
#endif
