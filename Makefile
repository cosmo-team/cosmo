# NOTE: needs boost, tclap, STXXL, and sdsl

CXX=g++#clang++-3.5#g++
CPP_FLAGS=-pipe -m64 -std=c++11 -pedantic-errors -W -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual \
					-Wunused -Wwrite-strings #\
					#-Werror
#-Wbool-conversions -Wshift-overflow -Wliteral-conversion \

DEP_PATH=/usr/local
INC=-I$(DEP_PATH)/include
LIB=-L$(DEP_PATH)/lib -L./
BOOST_FLAGS=-DBOOST_LOG_DYN_LINK -lboost_log -lboost_system -lboost_filesystem 
DEP_FLAGS=$(INC) $(LIB) $(BOOST_FLAGS)
DEBUG_FLAGS=-pg -gstabs
NDEBUG_FLAGS=-DNDEBUG
OPT_FLAGS=-O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing
NOPT_FLAGS=-O0

# Using Semantic Versioning: http://semver.org/
VERSION=0.5.1
CPP_FLAGS+=-DVERSION=\"$(VERSION)\"

k?=32
CPP_FLAGS+=-DK_LEN=$(k)

ifeq ($(optimise),0)
CPP_FLAGS+=$(NOPT_FLAGS)
else
CPP_FLAGS+=$(OPT_FLAGS)
endif

ifeq ($(debug),1)
CPP_FLAGS+=$(DEBUG_FLAGS)
else
CPP_FLAGS+=$(NDEBUG_FLAGS)
endif

ifeq ($(verbose),1)
CPP_FLAGS+=-DVERBOSE
endif

ifneq ($(revcomps),0)
CPP_FLAGS+=-DADD_REVCOMPS
endif

ifneq ($(dummies),0)
CPP_FLAGS+=-DALL_DUMMIES
endif

ifeq ($(varord),1)
CPP_FLAGS+=-DVAR_ORDER
endif

BUILD_REQS=debruijn_graph.hpp io.hpp io.o debug.hpp
ASSEM_REQS=debruijn_graph.hpp algorithm.hpp utility.hpp kmer.hpp
PACK_REQS=lut.hpp debug.hpp io.hpp io.o sort.hpp kmer.hpp dummies.hpp
BINARIES=cosmo-pack cosmo-build cosmo-benchmark # cosmo-test # cosmo-assemble

default: all

lut.hpp: make_lut.py
		python make_lut.py > lut.hpp

io.o: io.hpp io.cpp debug.hpp dummies.hpp kmer.hpp
		$(CXX) $(CPP_FLAGS) -c io.cpp

# TODO: Roll these all into one... "cosmo". Like git started off as multiple programs.
cosmo-pack: cosmo-pack.cpp $(PACK_REQS)
		$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(DEP_FLAGS) \
		-lstxxl -fopenmp -lsdsl #-lhpthread

cosmo-build: cosmo-build.cpp $(BUILD_REQS)
		$(CXX) $(CPP_FLAGS) -o $@ $< io.o $(DEP_FLAGS) -lsdsl

#cosmo-assemble: cosmo-assemble.cpp $(ASSEM_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
#		$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS) -lsdsl

cosmo-benchmark: cosmo-benchmark.cpp $(ASSEM_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
		$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS) -lsdsl

all: $(BINARIES)

catch.hpp:
	wget https://raw.githubusercontent.com/philsquared/Catch/master/single_include/catch.hpp

cosmo-test: cosmo-test.cpp catch.hpp $(wildcard *_test.cpp) $(wildcard $(subst _test.cpp,.hpp,$(wildcard *_test.cpp)))
	$(CXX) $(CPP_FLAGS) -o $@ $(filter-out %.hpp,$^) $(DEP_FLAGS) \
	-lstxxl -fopenmp -lsdsl

clean:
		rm -rf $(BINARIES) *.o *.dSYM
