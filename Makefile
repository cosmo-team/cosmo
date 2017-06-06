# NOTE: needs boost, tclap, STXXL, KMC, and sdsl

CXX=g++
CPP_FLAGS=-pipe -m64 -std=c++14 -pedantic-errors -W -Wall -Wextra -Wpointer-arith \
					-Wunused -Wwrite-strings #-Wcast-qual #-Werror

CLANG_WARNINGS=-Wbool-conversions -Wshift-overflow -Wliteral-conversion # CLANG ONLY
BOOST_PATH=3rd_party_inst/boost

DEP_PATH=3rd_party_inst
KMC_PATH=./3rd_party_src/KMC
INC=-isystem $(DEP_PATH)/include -isystem $(BOOST_PATH)/include
LIB=-L$(DEP_PATH)/lib -L./ -L$(BOOST_PATH)/lib
BOOST_FLAGS= -lboost_system -lboost_filesystem
DEP_FLAGS=$(INC) $(LIB) $(BOOST_FLAGS) -isystem $(KMC_PATH) -lsdsl -fopenmp #openmp is needed for logging
DEBUG_FLAGS=-pg -gstabs
NDEBUG_FLAGS=-DNDEBUG
OPT_FLAGS=-O3 -mmmx -msse -msse2 -msse3 -msse4 -msse4.2 -march=native -fno-strict-aliasing
NOPT_FLAGS=-O0

# Using Semantic Versioning: http://semver.org/
VERSION=0.7.0
BANNER='Copyright Alex Bowe (c) 2016'
CPP_FLAGS+=-DVERSION=\"$(VERSION)\" -DBANNER=\"$(BANNER)\"

k?=32
# TODO: quantize k
CPP_FLAGS+=-DK_LEN=$(k)

colors?=128
CPP_FLAGS+=-DNUM_COLS=$(colors)

ifeq ($(asm),1)
CPP_FLAGS+=-S -fverbose-asm
endif

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

KMC_OBJS=$(KMC_PATH)/kmc_api/kmc_file.o $(KMC_PATH)/kmc_api/kmer_api.o $(KMC_PATH)/kmc_api/mmer.o
BUILD_REQS=lut.hpp debug.hpp utility.hpp io.hpp sort.hpp kmer.hpp dummies.hpp debruijn_graph.hpp debruijn_graph_shifted.hpp pack-color.hpp cosmo-color-pd.hpp 
COLOR_REQS=colored_debruijn_graph.hpp io.hpp debug.hpp
BINARIES=cosmo-build cosmo-color cosmo-test cosmo-benchmark cosmo-benchmark-varord pack-color cosmo-color-pd cosmo-read-color transpose cosmo-merge

default: all

# Stores last compiled options to know whether we need to recompile when comp variables change
.PHONY: force
compiler_flags: force
	@echo $(CPP_FLAGS) | cmp -s - $@ || echo $(CPP_FLAGS) > $@

lut.hpp: make_lut.py
	python make_lut.py > lut.hpp

cosmo-build: cosmo-build.cpp $(BUILD_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS) -lstxxl

cosmo-color: cosmo-color.cpp $(BUILD_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

cosmo-merge: cosmo-merge.cpp $(BUILD_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)


cosmo-color-pd: cosmo-color-pd.cpp $(BUILD_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

cosmo-read-color: cosmo-read-color.cpp $(BUILD_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

transpose: transpose.cpp $(BUILD_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

pack-color: pack-color.cpp $(BUILD_REQS) compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(KMC_OBJS) $(DEP_FLAGS)

cosmo-benchmark: cosmo-benchmark.cpp $(BUILD_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
	$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS)

cosmo-benchmark-varord: cosmo-benchmark.cpp $(BUILD_REQS) wt_algorithm.hpp debruijn_hypergraph.hpp
	$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS) -DVAR_ORDER

#cosmo-test: cosmo-test.cpp catch.hpp $(wildcard *_test.cpp) $(wildcard $(subst _test.cpp,.hpp,$(wildcard *_test.cpp)))
#	$(CXX) $(CPP_FLAGS) -o $@ $(filter-out %.hpp,$^) $(DEP_FLAGS) -lstxxl -fopenmp -lsdsl

cosmo-test: debruijn_graph_test.cpp $(BUILD_REQS) bgl_sdb_adapter.hpp multi_bit_vector.hpp compiler_flags
	$(CXX) $(CPP_FLAGS) -o $@ $< $(DEP_FLAGS) -lboost_unit_test_framework

test: cosmo-test
	./cosmo-test

all: $(BINARIES)

clean:
		rm -rf $(BINARIES) *.o *.dSYM
