
#include "OutputWriter.hh"

#include "FilledObservable.hh"
#include "DataStructure.hh"
#include "JetrateDataStructure.hh"
#include "DifferentialDataStructure.hh"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include <map>
using std::map;
#include <iostream>
using std::cout;
using std::endl;

OutputWriter::OutputWriter( const string& filename ) {
  cout << "OutputWriter::OutputWriter: opening for writing: " << filename << endl;
  outputfile= new TFile( filename.c_str(), "RECREATE" );
}

OutputWriter::~OutputWriter() {
  cout << "OutputWriter::~OutputWriter: closing file" << endl;
  outputfile->Close();
  delete outputfile;
}

void OutputWriter::write( const vector<FilledObservable*>& vobs ) {
  for( size_t iobs= 0; iobs < vobs.size(); iobs++ ) {
    map<string,DataStructure*> data= vobs[iobs]->getData();
    for( map<string,DataStructure*>::iterator iter= data.begin();
	 iter != data.end(); iter++ ) {
      string txt= vobs[iobs]->getName() + " " + iter->first;
      DataStructure* ds= iter->second;
      DifferentialDataStructure* dds= dynamic_cast<DifferentialDataStructure*>( ds );
      JetrateDataStructure* jrds= dynamic_cast<JetrateDataStructure*>( ds );
      if( jrds ) {
	vector<Double_t> points= jrds->getPoints();
	vector<Double_t> values= jrds->getValues();
	vector<Double_t> errors= jrds->getErrors();
	Int_t n= points.size();
	Double_t* xerrors= new Double_t[n];
	for( Int_t i= 0; i < n; i++ ) {
	  xerrors[i]= 0.0;
	}
	TGraphErrors tge( n, &(points[0]), &(values[0]), xerrors, &(errors[0]) );
	tge.SetTitle( txt.c_str() );
	tge.SetMarkerStyle( 20 );
	tge.SetMarkerSize( 1.0 );
	cout << "OutputWriter::write: writing TGraphErrors " << txt << endl;
	tge.Write( txt.c_str() );
      } 
      else if( dds ) {
	vector<Double_t> binedges= dds->getBinedges();
	vector<Double_t> values= dds->getValues();
	vector<Double_t> errors= dds->getErrors();
	TH1D hist( txt.c_str(), txt.c_str(), binedges.size()-1, &(binedges[0]) );
	hist.SetContent( &(values[0]) );
	hist.SetError( &(errors[0]) );
	cout << "OutputWriter::write: writing TH1D " << txt << endl;
	hist.Write();
      }
      else {
	cout << "OutputWriter::write: dynamic_cast failed, no output" << endl;
      }
    }
  }
  return;
}

