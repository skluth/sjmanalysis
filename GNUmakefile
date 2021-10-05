
# Makefile for sjmanalysis package

CXX      = g++
LD       = $(CXX)
RC       = rootcint
OPT      = -g
CXXSTD   = -std=c++11
CXXFLAGS = -Wall -fPIC $(OPT) $(CXXSTD) -Wno-deprecated-declarations
FC       = gfortran
FFLAGS   = -fPIC $(OPT) 

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
# FASTJETLIBS = $(shell $(FASTJETCONFIG) --libs --plugins )
FASTJETLIBS = $(shell $(FASTJETCONFIG) --libs --plugins ) -lfastjetcontribfragile

ROOTCONFIG = $(HOME)/Downloads/root/root_v6.22.06/bin/root-config
ROOTINC = $(shell $(ROOTCONFIG) --noauxcflags --cflags )
ROOTLIBS = $(shell $(ROOTCONFIG) --libs ) -ltbb
ROOTLIBDIR = $(shell $(ROOTCONFIG) --libdir )

HEPMC2PATH = $(HOME)/qcd/hepmc/hepmc2.06.09/install
HEPMC2INC = -I$(HEPMC2PATH)/include
HEPMC2LIBS = -L$(HEPMC2PATH)/lib -lHepMC
HEPMC2LIBDIR = $(HEPMC2PATH)/lib

CPPFLAGS = $(ROOTINC) $(FASTJETINC) $(HEPMC2INC)

SRCS = LEPNtupleReader.cc TFastJet.cc Analysis.cc DataStructure.cc \
JetrateDataStructure.cc DifferentialDataStructure.cc MatrixDataStructure.cc \
Observable.cc ObsDifferential.cc ObsJetrate.cc ObsFastJetDiff.cc \
ObsPartonShower.cc ObsEEC.cc ObsGroomed.cc ObservableFactory.cc \
FilledObservable.cc Unfolder.cc BbbUnfolder.cc MtxUnfolder.cc OutputWriter.cc \
LEPThrustCalculator.cc LEPYnmCalculator.cc PxThrustCalculator.cc \
FastJetYcutCalculator.cc FastJetEminCalculator.cc FastJetRCalculator.cc \
FastJetPxConeRCalculator.cc FastJetPxConeEminCalculator.cc \
LEPYcutCalculator.cc AnalysisProcessor.cc SjmConfigParser.cc \
LEP1NtupleReader.cc LEP2NtupleReader.cc NtupleReader.cc \
HepMCRootReader.cc

# HMC2SRCS = runhepmc2.cc


# Fortran stuff for thrust
FSRCS = pxlth4.f


LIB = libNtupleReader.so

DICT = AnalysisDict.cc
DICTLIB = lib$(DICT:.cc=.so)
DICTSRCS = Analysis.cc TH1DAnalysisObject.cc TGEAnalysisObject.cc

DEPS = $(SRCS:.cc=.d) $(filter-out $(SRCS:.cc=.d), $(DICTSRCS:.cc=.d) )


all: testsjmanalysis runjob

# Compile Fortran
%.o : %.f
	$(FC) $(FFLAGS) -c -o $@ $<


$(DEPS): %.d: %.cc
	$(CXX) $(CPPFLAGS) $(CXXSTD) -MM $< -MF $@
-include $(DEPS)


$(LIB): $(SRCS:.cc=.o) $(FSRCS:.f=.o)
	$(CXX) -shared -Wl,--no-as-needed $(ROOTLIBS) $(FASTJETLIBS) $(HEPMC2LIBS) -o $@ $^


testsjmanalysis: testsjmanalysis.cc $(LIB)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(GINCS) -o $@ $^ $(GLIBS) $(ROOTLIBS) $(HEPMC2LIBS) -lboost_program_options
	LD_LIBRARY_PATH=$(PWD):$(ROOTLIBDIR):$(HEPMC2LIBDIR) ./$@



GenEventDataDict.cc: GenEventData.hh GenEventDataLinkDef.hh
	rootcint -f $@ -c $^

testHepMCRootReader: testHepMCRootReader.cc HepMCRootReader.cc GenEventDataDict.cc
	$(CXX) $(CXXFLAGS) $(GINCS) $(ROOTINC) -o $@ $^ $(GLIBS) $(ROOTLIBS)
	LD_LIBRARY_PATH=$(ROOTLIBDIR) ./$@


runjob: runjob.cc $(LIB)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(ROOTLIBS) $(HEPMC2LIBS) -lboost_program_options


runhepmc2: runhepmc2.cc $(LIB) GenEventDataDict.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(ROOTLIBS) $(HEPMC2LIBS) -lboost_program_options

testHepMC2Reader: testHepMC2Reader.cc HepMC2Reader.cc
	$(CXX) $(CXXFLAGS) $(GINCS) $(HEPMC2INC)  $(ROOTINC) -o $@ $^ $(GLIBS) $(HEPMC2LIBS) $(ROOTLIBS)
	LD_LIBRARY_PATH=$(ROOTLIBDIR):$(HEPMC2LIBDIR) ./$@

$(DICT): $(DICTSRCS:.cc=.hh) $(DICT:Dict.cc=LinkDef.h)
	$(RC) -f $@ -c $^

$(DICTLIB): $(DICT:.cc=.o) $(DICTSRCS:.cc=.o)
	$(CXX) -shared -Wl,--no-as-needed $(ROOTLIBS) $(FASTJETLIBS) -o $@ $^


clean:
	rm -f $(SRCS:.cc=.o) $(FSRCS:.f=.o) $(LIB) $(DEPS) testsjmanalysis testsjmanalysis.o runjob $(DICT) $(DICTLIB) $(DICTSRCS:.cc=.o)

