ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

CXX       := g++
CXXFLAGS  += -std=c++17 -Wall -Wshadow -Warray-bounds -Wmissing-field-initializers -fPIC $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS) -O3
LDFLAGS   += -std=c++17 $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS) -O3

# These are the compiled files from the GenieAnalysis "library" (see later) needed for apps
GA_LIB_OBJECTS := GenieAnalysis/GenieAnalysis.o \
				  GenieAnalysis/GAAutoHistograms.o \
				  GenieAnalysis/GACLAS6Common.o \
				  GenieAnalysis/GACLAS6MC.o \
				  GenieAnalysis/GACLAS6Data.o \
				  GenieAnalysis/GACLAS6FinalFST.o \
				  GenieAnalysis/Fiducial.o \
				  GenieAnalysis/misc.o


# Executables
default: bin bin/genie_analysis_demo

# Folder where to put all compiled executables
bin:
	mkdir bin

bin/genie_analysis_demo: bin apps/genie_analysis_demo.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $(filter-out bin, $^) -I . $(LDFLAGS)

# Apps from Jan Kocka's SHP done in spring 2023
bin/final_truth_FST: bin apps/jkocka_spring2023/final_truth_FST.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $(filter-out bin, $^) -I . $(LDFLAGS)

bin/misc: bin apps/jkocka_spring2023/misc.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $(filter-out bin, $^) -I . $(LDFLAGS)

bin/genie_analysis_data: bin apps/jkocka_spring2023/genie_analysis_data.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $(filter-out bin, $^) -I . $(LDFLAGS)

bin/genie_analysis_1p: bin apps/jkocka_spring2023/genie_analysis_1p.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $(filter-out bin, $^) -I . $(LDFLAGS)

bin/genie_analysis_FS_transparency: bin apps/jkocka_spring2023/genie_analysis_FS_transparency.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $(filter-out bin, $^) -I . $(LDFLAGS)

bin/genie_analysis_1pi1nuc: bin apps/jkocka_spring2023/genie_analysis_1pi1nuc.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $(filter-out bin, $^) -I . $(LDFLAGS)

bin/genie_analysis_deltas: bin apps/jkocka_spring2023/genie_analysis_deltas.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $(filter-out bin, $^) -I . $(LDFLAGS)

bin/genie_analysis_lucass_cuts: bin apps/jkocka_spring2023/genie_analysis_lucass_cuts.cpp $(GA_LIB_OBJECTS)
	$(CXX) -o $@ $(filter-out bin, $^) -I . $(LDFLAGS)

# Groups of eecutables - not really needed, but nice to organize stuff
jkocka_spring2023: bin/genie_analysis_deltas bin/genie_analysis_1pi1nuc bin/genie_analysis_1p \
	bin/genie_analysis_FS_transparency bin/genie_analysis_data bin/final_truth_FST

all: bin bin/genie_analysis_demo jkocka_spring2023


# GenieAnalysis "library"
#
# Fiducial is taken from original and not touched as long as it works! (so disable the ugly warnings, remove if modified)
GenieAnalysis/Fiducial.o: GenieAnalysis/Fiducial.cpp GenieAnalysis/Fiducial.h
	$(CXX) -o $@ -c $< $(CXXFLAGS) -w
# everything else
GenieAnalysis/%.o: GenieAnalysis/%.cpp GenieAnalysis/%.h
	$(CXX) -o $@ -c $< $(CXXFLAGS)


clean:
	@rm -rf GenieAnalysis/*.o bin

debug: CXXFLAGS += -pg
debug: LDFLAGS += -pg
debug: all
