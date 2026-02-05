
# Makefile for sjmanalysis package

CXX      = g++
LD       = $(CXX)
RC       = rootcint
OPT      = -g
CXXSTD   = -std=c++17
CXXFLAGS = -Wall -fPIC $(OPT) $(CXXSTD) -Wno-deprecated-declarations
FC       = gfortran
FFLAGS   = -fPIC $(OPT)

# Googletest
GTESTPATH = /opt/googletest/googletest-1.16.0
GINCS = -I $(GTESTPATH)/googletest/include -I $(GTESTPATH)/googlemock/include
GLIBS = -L $(GTESTPATH)/lib -lgmock -lgtest -lpthread

# Fastjet
FASTJETCONFIG = /opt/fastjet/fastjet-3.5.1/bin/fastjet-config
FASTJETINC = $(shell $(FASTJETCONFIG) --cxxflags )
FASTJETPATH = $(shell $(FASTJETCONFIG) --prefix )
FASTJETLIBDIR = $(FASTJETPATH)/lib
FASTJETLIBS = -L$(FASTJETLIBDIR) -lfastjetplugins -lsiscone_spherical -lsiscone -lRecursiveTools -lfastjettools -lfastjet -lgfortran -lm -lquadmath

# ROOT
ROOTCONFIG = $(shell which root-config )
ROOTINC = $(shell $(ROOTCONFIG) --noauxcflags --cflags )
# Recent ROOT >= 6.24/06 root-config emits libraries inconsistent :(
BADLIBS = -lROOTDataFrame -lROOTVecOps
ROOTLIBS = $(filter-out $(BADLIBS), $(shell $(ROOTCONFIG) --libs ))
ROOTLIBDIR = $(shell $(ROOTCONFIG) --libdir )

CPPFLAGS = $(ROOTINC) $(FASTJETINC)

SRCS = LEPNtupleReader.cc TFastJet.cc Analysis.cc DataStructure.cc \
JetrateDataStructure.cc DifferentialDataStructure.cc MatrixDataStructure.cc \
Observable.cc ObsDifferential.cc ObsJetrate.cc ObsFastJetDiff.cc \
ObsPartonShower.cc ObsEEC.cc ObsJEEC.cc ObsDEEC.cc ObsGroomed.cc ObservableFactory.cc \
FilledObservable.cc Unfolder.cc BbbUnfolder.cc MtxUnfolder.cc OutputWriter.cc \
LEPThrustCalculator.cc LEPYnmCalculator.cc PxThrustCalculator.cc \
FastJetYcutCalculator.cc FastJetEminCalculator.cc FastJetRCalculator.cc \
FastJetPxConeRCalculator.cc FastJetPxConeEminCalculator.cc \
LEPYcutCalculator.cc AnalysisProcessor.cc SjmConfigParser.cc \
LEP1NtupleReader.cc LEP2NtupleReader.cc NtupleReader.cc \
HepMCRootReader.cc

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

# Dependencies
$(DEPS): %.d: %.cc
	$(CXX) $(CPPFLAGS) $(CXXSTD) -MM $< -MF $@
-include $(DEPS)

$(LIB): $(SRCS:.cc=.o) $(FSRCS:.f=.o)
	$(CXX) -shared -Wl,--no-as-needed $(ROOTLIBS) $(FASTJETLIBS) $(HEPMC2LIBS) -o $@ $^

testsjmanalysis: testsjmanalysis.cc $(LIB)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(GINCS) -o $@ $^ $(GLIBS) $(ROOTLIBS) $(FASTJETLIBS) -lboost_program_options
	LD_LIBRARY_PATH=$(PWD):$(ROOTLIBDIR):$(FASTJETLIBDIR) ./$@

GenEventDataDict.cc: GenEventData.hh GenEventDataLinkDef.hh
	rootcint -f $@ -c $^

testHepMCRootReader: testHepMCRootReader.cc HepMCRootReader.cc GenEventDataDict.cc $(LIB)
	$(CXX) $(CXXFLAGS) $(GINCS) $(ROOTINC) -o $@ $^ $(GLIBS) $(ROOTLIBS) $(FASTJETLIBS) -lboost_program_options
	LD_LIBRARY_PATH=$(PWD):$(ROOTLIBDIR):$(FASTJETLIBDIR) ./$@

runjob: runjob.cc $(LIB)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(ROOTLIBS) $(FASTJETLIBS) -lboost_program_options

runhepmc2: runhepmc2.cc $(LIB) GenEventDataDict.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(ROOTLIBS) $(FASTJETLIBS) -lboost_program_options


$(DICT): $(DICTSRCS:.cc=.hh) $(DICT:Dict.cc=LinkDef.h)
	$(RC) -f $@ -c $^
$(DICTLIB): $(DICT:.cc=.o) $(DICTSRCS:.cc=.o)
	$(CXX) -shared -Wl,--no-as-needed $(ROOTLIBS) $(FASTJETLIBS) -o $@ $^


# Create cpython binding for YKERN for tests, use "import ylcus" in python
yclus.cpython-38-x86_64-linux-gnu.so: yclus.f
	f2py3 -c yclus.f -m yclus

clean:
	rm -f $(SRCS:.cc=.o) $(FSRCS:.f=.o) $(LIB) $(DEPS) testsjmanalysis testsjmanalysis.o runjob $(DICT) $(DICTLIB) $(DICTSRCS:.cc=.o)

