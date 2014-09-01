# NOTE: needs boost, tclap, and sdsl

INCLUDE_PATH=/usr/local/include
LIBRARY_PATH=/usr/local/lib
CXX=g++
CPP_FLAGS=-m64 -std=c++0x -pedantic-errors -W -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual \
					-Wunused -Wstrict-prototypes -Wmissing-prototypes -Wwrite-strings -Werror
DEP_FLAGS=-I/usr/local/include -L/usr/local/lib -lsdsl
DEBUG_FLAGS=-g -O0
# The MMX and SSE flags can be safely disabled
RELEASE_FLAGS=-O3 -DNDEBUG -mmmx -msse -msse2 -msse3 -msse4 -march=native

# DEFINITIONS
# Using Semantic Versioning: http://semver.org/
VERSION=0.4.2
CPP_FLAGS+=-DVERSION=\"$(VERSION)\"

ifeq ($(debug),1)
CPP_FLAGS+=$(DEBUG_FLAGS)
else
CPP_FLAGS+=$(RELEASE_FLAGS)
endif

ifneq ($(revcomps),0)
CPP_FLAGS+=-DADD_REVCOMPS
endif

ifneq ($(dummies),0)
CPP_FLAGS+=-DALL_DUMMIES
endif

DBG_REQS=debruijn_graph.hpp utility.hpp io.hpp io.o debug.h
PACK_REQS=lut.hpp debug.h io.hpp io.o sort.hpp kmer.hpp dummies.hpp
BINARIES=pack-edges cosmo-build cosmo-assemble

default: all

lut.hpp: make_lut.py
		python make_lut.py > lut.hpp

io.o: io.hpp io.cpp debug.h dummies.hpp kmer.hpp
		$(CXX) $(CPP_FLAGS) -c io.cpp

# TODO: Roll these all into one... "cosmo"
pack-edges: pack-edges.cpp $(PACK_REQS)
		$(CXX) $(CPP_FLAGS) -o $@ $< io.o

cosmo-build: cosmo-build.cpp $(DBG_REQS)
		$(CXX) $(CPP_FLAGS) $(DEP_FLAGS) -o $@ $< io.o

cosmo-assemble: cosmo-assemble.cpp $(DBG_REQS) algorithm.hpp
		$(CXX) $(CPP_FLAGS) $(DEP_FLAGS) -o $@ $<

all: $(BINARIES)

clean:
		rm -rf $(BINARIES) *.o *.dSYM
