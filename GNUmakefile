
# Makefile for sjmanalysis package

CXX      = g++
LD       = $(CXX)
RC       = rootcint
OPT      = -g
CXXFLAGS = -Wall -fPIC $(OPT) -std=c++11

#GTESTPATH = $(HOME)/Downloads/googletest/googletest-master/googletest
#GMOCKPATH = $(HOME)/Downloads/googletest/googletest-master/googlemock
# GMOCKPATH = $(HOME)/Downloads/googletest/googlemock
#GINCS = -I $(GMOCKPATH)/include -I $(GTESTPATH)/include
#GLIBS = -L $(GMOCKPATH) -l gmock -L $(GTESTPATH) -lgtest -lpthread
GTESTPATH = $(HOME)/Downloads/googletest
GINCS = -I $(GTESTPATH)/include
GLIBS = -L $(GTESTPATH)/lib -l gmock -l gtest -l pthread

FASTJETCONFIG = $(HOME)/qcd/fastjet/fastjet-3.1.3/install/bin/fastjet-config
FASTJETINC = $(shell $(FASTJETCONFIG) --cxxflags )
FASTJETLIBS = $(shell $(FASTJETCONFIG) --libs --plugins )

ROOTINC = $(shell root-config --noauxcflags --cflags )
ROOTLIBS = $(shell root-config --libs )
ROOTLIBDIR = $(shell root-config --libdir )

CPPFLAGS = $(ROOTINC) $(FASTJETINC)

SRCS = NtupleReader.cc TFastJet.cc Analysis.cc \
JetrateDataStructure.cc DifferentialDataStructure.cc MatrixDataStructure.cc \
Observable.cc ObsDifferential.cc ObsJetrate.cc ObsFastJetDiff.cc \
ObsPartonShower.cc ObservableFactory.cc \
FilledObservable.cc Unfolder.cc OutputWriter.cc \
ThrustCalculator.cc YnmdCalculator.cc YnmjCalculator.cc \
FastJetYcutCalculator.cc FastJetEminCalculator.cc FastJetRCalculator.cc \
FastJetPxConeRCalculator.cc FastJetPxConeEminCalculator.cc \
YcutCalculator.cc AnalysisProcessor.cc SjmConfigParser.cc \
LEP1NtupleReader.cc LEP2NtupleReader.cc

LIB = libNtupleReader.so
DEPS = $(SRCS:.cc=.d)

all: testsjmanalysis runjob

$(DEPS): %.d: %.cc
	$(CXX) $(CPPFLAGS) -MM $< -MF $@
-include $(DEPS)

$(LIB): $(SRCS:.cc=.o)
	$(CXX) -shared -Wl,--no-as-needed $(ROOTLIBS) $(FASTJETLIBS) -o $@ $^

testsjmanalysis: testsjmanalysis.cc $(LIB)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(GINCS) -o $@ $^ $(GLIBS) $(ROOTLIBS) -lboost_program_options
	LD_LIBRARY_PATH=$(PWD):$(ROOTLIBDIR) ./$@

runjob: runjob.cc $(LIB)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(ROOTLIBS) -lboost_program_options

clean:
	rm -f $(SRCS:.cc=.o) $(LIB) $(DEPS) testsjmanalysis testsjmanalysis.o runjob
