
# Makefile for sjmanalysis package

CXX      = g++
LD       = $(CXX)
RC       = rootcint
OPT      = -g
CXXSTD   = -std=c++11
CXXFLAGS = -Wall -fPIC $(OPT) $(CXXSTD) -Wno-deprecated-declarations

#GTESTPATH = $(HOME)/Downloads/googletest/googletest-master/googletest
#GMOCKPATH = $(HOME)/Downloads/googletest/googletest-master/googlemock
# GMOCKPATH = $(HOME)/Downloads/googletest/googlemock
#GINCS = -I $(GMOCKPATH)/include -I $(GTESTPATH)/include
#GLIBS = -L $(GMOCKPATH) -l gmock -L $(GTESTPATH) -lgtest -lpthread
GTESTPATH = $(HOME)/Downloads/googletest
GINCS = -I $(GTESTPATH)/include
GLIBS = -L $(GTESTPATH)/lib -l gmock -l gtest -l pthread

FASTJETCONFIG = $(HOME)/qcd/fastjet/fastjet-3.3.0/install/bin/fastjet-config
FASTJETINC = $(shell $(FASTJETCONFIG) --cxxflags )
FASTJETLIBS = $(shell $(FASTJETCONFIG) --libs --plugins )

ROOTCONFIG = $(HOME)/Downloads/root/root_v6.12.04/bin/root-config
ROOTINC = $(shell $(ROOTCONFIG) --noauxcflags --cflags )
ROOTLIBS = $(shell $(ROOTCONFIG) --libs )
ROOTLIBDIR = $(shell $(ROOTCONFIG) --libdir )

CPPFLAGS = $(ROOTINC) $(FASTJETINC)

SRCS = NtupleReader.cc TFastJet.cc Analysis.cc DataStructure.cc \
JetrateDataStructure.cc DifferentialDataStructure.cc MatrixDataStructure.cc \
Observable.cc ObsDifferential.cc ObsJetrate.cc ObsFastJetDiff.cc \
ObsPartonShower.cc ObsEEC.cc ObservableFactory.cc \
FilledObservable.cc Unfolder.cc BbbUnfolder.cc MtxUnfolder.cc OutputWriter.cc \
ThrustCalculator.cc YnmCalculator.cc \
FastJetYcutCalculator.cc FastJetEminCalculator.cc FastJetRCalculator.cc \
FastJetPxConeRCalculator.cc FastJetPxConeEminCalculator.cc \
YcutCalculator.cc AnalysisProcessor.cc SjmConfigParser.cc \
LEP1NtupleReader.cc LEP2NtupleReader.cc

LIB = libNtupleReader.so

DICT = AnalysisDict.cc
DICTLIB = lib$(DICT:.cc=.so)
DICTSRCS = Analysis.cc TH1DAnalysisObject.cc TGEAnalysisObject.cc

DEPS = $(SRCS:.cc=.d) $(filter-out $(SRCS:.cc=.d), $(DICTSRCS:.cc=.d) )

all: testsjmanalysis runjob

$(DEPS): %.d: %.cc
	$(CXX) $(CPPFLAGS) $(CXXSTD) -MM $< -MF $@
-include $(DEPS)

$(LIB): $(SRCS:.cc=.o)
	$(CXX) -shared -Wl,--no-as-needed $(ROOTLIBS) $(FASTJETLIBS) -o $@ $^

testsjmanalysis: testsjmanalysis.cc $(LIB)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(GINCS) -o $@ $^ $(GLIBS) $(ROOTLIBS) -lboost_program_options
	LD_LIBRARY_PATH=$(PWD):$(ROOTLIBDIR) ./$@

runjob: runjob.cc $(LIB)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(ROOTLIBS) -lboost_program_options

$(DICT): $(DICTSRCS:.cc=.hh) $(DICT:Dict.cc=LinkDef.h)
	$(RC) -f $@ -c $^

$(DICTLIB): $(DICT:.cc=.o) $(DICTSRCS:.cc=.o)
	$(CXX) -shared -Wl,--no-as-needed $(ROOTLIBS) $(FASTJETLIBS) -o $@ $^


clean:
	rm -f $(SRCS:.cc=.o) $(LIB) $(DEPS) testsjmanalysis testsjmanalysis.o runjob $(DICT) $(DICTLIB) $(DICTSRCS:.cc=.o)

