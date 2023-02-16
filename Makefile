ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

CXX       := g++
CXXFLAGS  += -std=c++17 -Wall -Wshadow -Warray-bounds -Wmissing-field-initializers -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)


# Executables
genie_analysis_lucass_cuts: GenieAnalysis/GenieAnalysis.o GenieAnalysis/GenieAnalysisAuto.o GenieAnalysis/misc.h GenieAnalysis/Fiducial.o genie_analysis_lucass_cuts.cpp
	$(CXX) -o $@ $^ $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)

genie_analysis_demo: GenieAnalysis/GenieAnalysis.o GenieAnalysis/GenieAnalysisAuto.o GenieAnalysis/misc.h genie_analysis_demo.cpp
	$(CXX) -o $@ $^ $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)


# GenieAnalysis "library"
GenieAnalysis/%.o: GenieAnalysis/%.cpp GenieAnalysis/%.h
	$(CXX) -o $@ -c $< -O2 $(CXXFLAGS) $(INCLUDES)


clean:
	@rm -rf GenieAnalysis/*.o *.o genie_analysis_demo
