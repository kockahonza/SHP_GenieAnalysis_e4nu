ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

CXX       := g++
CXXFLAGS  += -std=c++17 -Wall -Wshadow -Warray-bounds -Wmissing-field-initializers -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

%.o: %.C %.h
	$(CXX) -c $< -O2 $(CXXFLAGS) $(INCLUDES)

genie_analysis_demo: genie_analysis.o
	$(CXX) -o genie_analysis_demo genie_analysis_demo.C $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS)

clean:
	@rm -rf *.o genie_analysis_demo
