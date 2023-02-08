ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

CXX       := g++
CXXFLAGS  += -std=c++17 -Wall -Wshadow -Warray-bounds -Wmissing-field-initializers -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

GenieAnalysis/%.o: GenieAnalysis/%.cpp GenieAnalysis/%.h
	$(CXX) -o $@ -c $< -O2 $(CXXFLAGS) $(INCLUDES)

genie_analysis_demo: GenieAnalysis/GenieAnalysis.o GenieAnalysis/misc.h genie_analysis_demo.cpp
	$(CXX) -o genie_analysis_demo genie_analysis_demo.cpp GenieAnalysis/GenieAnalysis.o $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)

clean:
	@rm -rf GenieAnalysis/*.o *.o genie_analysis_demo
