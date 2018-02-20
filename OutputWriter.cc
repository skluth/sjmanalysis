
#include "OutputWriter.hh"

#include "FilledObservable.hh"
#include "DataStructure.hh"
#include "JetrateDataStructure.hh"
#include "DifferentialDataStructure.hh"
#include "MatrixDataStructure.hh"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include <map>
using std::map;
#include <iostream>
using std::cout;
using std::endl;
#include <stdexcept>
using std::logic_error;

OutputWriter::OutputWriter( const string& filename ) {
  cout << "OutputWriter::OutputWriter: opening for writing: " << filename << endl;
  outputfile= new TFile( filename.c_str(), "RECREATE" );
  if( not outputfile->IsOpen() ) {
    string txt= "OutputWriter::OutputWriter: file not open: " + filename;
    throw std::logic_error( txt );
  }
}

OutputWriter::~OutputWriter() {
  cout << "OutputWriter::~OutputWriter: closing file" << endl;
  outputfile->Close();
  delete outputfile;
}

void
OutputWriter::writeJetrate( const JetrateDataStructure* jrds,
			    const string& txt ) {
  vector<Double_t> points= jrds->getPoints();
  vector<Double_t> values= jrds->getValues();
  vector<Double_t> errors= jrds->getErrors();
  Int_t n= points.size();
  Double_t xerrors[n];
  for( Int_t i= 0; i < n; i++ ) xerrors[i]= 0.0;
  TGraphErrors tge( n, &(points[0]), &(values[0]), xerrors, &(errors[0]) );
  if( not jrds->getNormalised() ) {
    tge.SetPoint( n, 99.0, jrds->getNEvents() );
    tge.SetPointError( n, 0.0, 0.0 );
  }
  tge.SetTitle( txt.c_str() );
  tge.SetMarkerStyle( 20 );
  tge.SetMarkerSize( 1.0 );
  cout << "OutputWriter::writeJetrate: writing TGraphErrors " << txt << endl;
  tge.Write( txt.c_str() );
  return;
}

void
OutputWriter::writeDifferentialDistribution( const DifferentialDataStructure* dds, 
					     const string& txt ) {
  vector<Double_t> binedges= dds->getBinedges();
  vector<Double_t> values= dds->getValues();
  vector<Double_t> errors= dds->getErrors();
  TH1D hist( txt.c_str(), txt.c_str(), binedges.size()-1, &(binedges[0]) );
  hist.SetContent( &(values[0]) );
  hist.SetError( &(errors[0]) );
  hist.SetEntries( dds->getNEvents() );
  cout << "OutputWriter::writeDifferentialDistribution: writing TH1D " << txt << endl;
  hist.Write();
  return;
}

void OutputWriter::writeMatrix( MatrixDataStructure* mds, 
				const string& txt ) {
  vector<Double_t> binedges= mds->getBinedges();
  Int_t nbin= binedges.size()-1;
  TH2D hist( txt.c_str(), txt.c_str(), nbin, &(binedges[0]), nbin, &(binedges[0]) );
  for( Int_t ibin= 0; ibin < nbin+2; ibin++ ) {
    for( Int_t jbin= 0; jbin < nbin+2; jbin++ ) {
      hist.SetBinContent( jbin, ibin, mds->getElement( jbin, ibin ) );
    }
  }
  hist.SetEntries( mds->getNEvents() );
  cout << "OutputWriter::writeMatrix: writing TH2D " << txt << endl;
  hist.Write();
  return;
}

void OutputWriter::write( const vector<FilledObservable*>& vobs ) {
  for( size_t iobs= 0; iobs < vobs.size(); iobs++ ) {
    map<string,DataStructure*> data= vobs[iobs]->getData();
    for( map<string,DataStructure*>::iterator iter= data.begin();
	 iter != data.end(); iter++ ) {
      string txt= vobs[iobs]->getName() + " " + iter->first;
      DataStructure* ds= iter->second;
      JetrateDataStructure* jrds= dynamic_cast<JetrateDataStructure*>( ds );
      DifferentialDataStructure* dds= dynamic_cast<DifferentialDataStructure*>( ds );
      if( jrds ) writeJetrate( jrds, txt );
      else if( dds ) writeDifferentialDistribution( dds, txt );
      else throw logic_error( "OutputWriter::write: wrong class" );
    }
    map<string,MatrixDataStructure*> matrices= vobs[iobs]->getMatrices();
    for( map<string,MatrixDataStructure*>::iterator iter= matrices.begin();
	 iter != matrices.end(); iter++ ) {
      MatrixDataStructure* mds= iter->second;
      string txt= vobs[iobs]->getName() + " " + iter->first;
      writeMatrix( mds, txt );
    }
  }
  return;
}

