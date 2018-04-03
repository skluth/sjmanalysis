
#include "OutputWriter.hh"

#include "FilledObservable.hh"
#include "DataStructure.hh"
#include "JetrateDataStructure.hh"
#include "DifferentialDataStructure.hh"
#include "MatrixDataStructure.hh"
#include "Normalisation.hh"

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
using std::runtime_error;

using std::string;
using std::vector;

OutputWriter::OutputWriter( const string & filename ) {
  cout << "OutputWriter::OutputWriter: opening for writing: " << filename << endl;
  outputfile= new TFile( filename.c_str(), "RECREATE" );
  if( not outputfile->IsOpen() ) {
    string txt= "OutputWriter::OutputWriter: file not open: " + filename;
    throw std::runtime_error( txt );
  }
}

OutputWriter::~OutputWriter() {
  cout << "OutputWriter::~OutputWriter: closing file" << endl;
  outputfile->Close();
  delete outputfile;
}

void
OutputWriter::writeJetrate( const JetrateDataStructure* jrds,
			    const string & txt ) {
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
  SetIsNormalised( tge, jrds->getNormalised() );
  cout << "OutputWriter::writeJetrate: writing TGraphErrors " << txt << endl;
  tge.Write( txt.c_str() );
  return;
}

void
OutputWriter::writeDifferentialDistribution( const DifferentialDataStructure* dds, 
					     const string & txt ) {
  vector<Double_t> binedges= dds->getBinedges();
  vector<Double_t> values= dds->getValues();
  vector<Double_t> errors= dds->getErrors();
  TH1D hist( txt.c_str(), txt.c_str(), binedges.size()-1, &(binedges[0]) );
  hist.SetContent( &(values[0]) );
  hist.SetError( &(errors[0]) );
  hist.SetEntries( dds->getNEvents() );
  SetIsNormalised( hist, dds->getNormalised() );
  cout << "OutputWriter::writeDifferentialDistribution: writing TH1D " << txt << endl;
  hist.Write();
  return;
}

void OutputWriter::writeMatrix( MatrixDataStructure* mds, 
				const string & txt,
				const string & prefix ) {
  vector<Double_t> binedges= mds->getBinedges();
  Int_t nbin= binedges.size()-1;
  string keytxt= prefix+" "+txt;
  TH2D hist( keytxt.c_str(), txt.c_str(), nbin, &(binedges[0]), nbin, &(binedges[0]) );
  for( Int_t ibin= 0; ibin < nbin+2; ibin++ ) {
    for( Int_t jbin= 0; jbin < nbin+2; jbin++ ) {
      hist.SetBinContent( jbin, ibin, mds->getElement( jbin, ibin ) );
    }
  }
  hist.SetEntries( mds->getNEvents() );
  cout << "OutputWriter::writeMatrix: writing TH2D " << keytxt << endl;
  hist.Write();
  return;
}

void OutputWriter::write( const vector<FilledObservable*> & vobs ) {
  for( const FilledObservable* obs : vobs ) {
    map<string,DataStructure*> data= obs->getData();
    for( const auto & keyValue : data ) {
      string txt= obs->getName() + " " + keyValue.first;
      DataStructure* ds= keyValue.second;
      JetrateDataStructure* jrds= dynamic_cast<JetrateDataStructure*>( ds );
      DifferentialDataStructure* dds= dynamic_cast<DifferentialDataStructure*>( ds );
      if( jrds ) writeJetrate( jrds, txt );
      else if( dds ) writeDifferentialDistribution( dds, txt );
      else throw logic_error( "OutputWriter::write: wrong class" );
    }
    map<string,MatrixDataStructure*> migrationMatrices= obs->getMigrationMatrices();
    for( const auto & keyValue : migrationMatrices ) {
      MatrixDataStructure* mds= keyValue.second;
      string txt= obs->getName() + " " + keyValue.first;
      writeMatrix( mds, txt, "migr" );
    }
    map<string,MatrixDataStructure*> errorMatrices= obs->getErrorMatrices();
    for( const auto & keyValue : errorMatrices ) {
      MatrixDataStructure* mds= keyValue.second;
      if( mds != 0 ) {
	string txt= obs->getName() + " " + keyValue.first;
	writeMatrix( mds, txt, "errm" );
      }
    }
  }
  return;
}
