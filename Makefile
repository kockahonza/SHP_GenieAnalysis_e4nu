ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

GA_LIB_OBJECTS := GenieAnalysis/GenieAnalysis.o \
				  GenieAnalysis/GAAutoHistograms.o \
				  GenieAnalysis/GACLAS6MC.o \
				  GenieAnalysis/GACLAS6Data.o \
				  GenieAnalysis/Fiducial.o \
				  GenieAnalysis/misc.o

CXX       := g++
CXXFLAGS  += -std=c++17 -Wall -Wshadow -Warray-bounds -Wmissing-field-initializers -fPIC $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS) -O3
LDFLAGS   += -std=c++17 $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS) -O3


# Executables
all: genie_analysis_deltas genie_analysis_1pi1nuc genie_analysis_FS_transparency genie_analysis_data

genie_analysis_data: genie_analysis_data.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

genie_analysis_1p: genie_analysis_1p.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

genie_analysis_FS_transparency: genie_analysis_FS_transparency.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

genie_analysis_1pi1nuc: genie_analysis_1pi1nuc.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

genie_analysis_deltas: genie_analysis_deltas.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

genie_analysis_lucass_cuts: genie_analysis_lucass_cuts.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)

genie_analysis_demo: genie_analysis_demo.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS)


# GenieAnalysis "library"
#
# Fiducial is taken from original and not touched as long as it works! (so disable the ugly warnings, remove if modified)
GenieAnalysis/Fiducial.o: GenieAnalysis/Fiducial.cpp GenieAnalysis/Fiducial.h
	$(CXX) -o $@ -c $< $(CXXFLAGS) -w
# everything else
GenieAnalysis/%.o: GenieAnalysis/%.cpp GenieAnalysis/%.h
	$(CXX) -o $@ -c $< $(CXXFLAGS)


clean:
	@rm -rf GenieAnalysis/*.o *.o \
		genie_analysis_demo \
		genie_analysis_lucass_cuts \
		genie_analysis_deltas \
		genie_analysis_1pi1nuc \
		genie_analysis_FS_transparency \
		genie_analysis_1p \
		genie_analysis_data

debug: CXXFLAGS += -pg
debug: LDFLAGS += -pg
debug: all
