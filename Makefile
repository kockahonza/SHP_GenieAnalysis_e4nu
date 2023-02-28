ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

GA_LIB_OBJECTS := GenieAnalysis/GenieAnalysis.o GenieAnalysis/misc.o GenieAnalysis/GenieAnalysisAuto.o GenieAnalysis/Fiducial.o GenieAnalysis/GenieAnalysisOriginalCuts.o

CXX       := g++
CXXFLAGS  += -std=c++17 -Wall -Wshadow -Warray-bounds -Wmissing-field-initializers -fPIC $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS) -O3
LDFLAGS   += -std=c++17 $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)


# Executables
genie_analysis_jan: genie_analysis_jan.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

genie_analysis_lucass_cuts: genie_analysis_lucass_cuts.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

genie_analysis_demo: genie_analysis_demo.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)


ALL: genie_analysis_demo genie_analysis_lucass_cuts genie_analysis_jan

# GenieAnalysis "library"
#
# Fiducial is taken from original and not touched as long as it works! (so disable the ugly warnings, remove if modified)
GenieAnalysis/Fiducial.o: GenieAnalysis/Fiducial.cpp GenieAnalysis/Fiducial.h
	$(CXX) -o $@ -c $< $(CXXFLAGS) -w
# everything else
GenieAnalysis/%.o: GenieAnalysis/%.cpp GenieAnalysis/%.h
	$(CXX) -o $@ -c $< $(CXXFLAGS)


clean:
	@rm -rf GenieAnalysis/*.o *.o genie_analysis_demo genie_analysis_lucass_cuts
