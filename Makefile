ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

GA_LIB_OBJECTS := GenieAnalysis/GenieAnalysis.o GenieAnalysis/misc.o GenieAnalysis/GenieAnalysisAuto.o GenieAnalysis/Fiducial.o GenieAnalysis/GenieAnalysisOriginalCuts.o

CXX       := g++
CXXFLAGS  += -std=c++17 -Wall -Wshadow -Warray-bounds -Wmissing-field-initializers -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)


# Executables
genie_analysis_lucass_cuts: genie_analysis_lucass_cuts.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $^ $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)

genie_analysis_demo: genie_analysis_demo.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $^ $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)


# GenieAnalysis "library"
GenieAnalysis/%.o: GenieAnalysis/%.cpp GenieAnalysis/%.h
	$(CXX) -o $@ -c $< -O2 $(CXXFLAGS) $(INCLUDES)


clean:
	@rm -rf GenieAnalysis/*.o *.o genie_analysis_demo genie_analysis_lucass_cuts
