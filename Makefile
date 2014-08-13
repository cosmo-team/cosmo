CXX=g++
CPP_FLAGS=-m64 -std=c++0x -pedantic-errors -W -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual \
					-Wunused -Wstrict-prototypes -Wmissing-prototypes -Wwrite-strings -Werror -march=native
DEBUG_FLAGS=-g -O0
RELEASE_FLAGS=-O3 -DNDEBUG -msse -msse2 -msse3 -msse4
REQS=convert_kmers.cpp lut.hpp debug.h nanotime.h io.o sort.hpp kmer.hpp dummies.hpp
COMPILE=$(CXX) $(CPP_FLAGS)

default: all

lut.hpp: make_lut.py
		python make_lut.py > lut.hpp

io.o: io.hpp io.cpp debug.h dummies.hpp kmer.hpp
		$(COMPILE) $(RELEASE_FLAGS) -c io.cpp

all: $(REQS)
		$(COMPILE) $(RELEASE_FLAGS) -o convert_kmers convert_kmers.cpp io.o

debug: $(REQS)
		$(COMPILE) $(DEBUG_FLAGS) -o convert_kmers convert_kmers.cpp io.o

clean:
		rm -rf convert_kmers *.o *.dSYM lut.hpp
