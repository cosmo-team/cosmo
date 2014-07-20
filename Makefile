CPP=g++
CPP_FLAGS=-m64 -std=c++0x -pedantic-errors -W -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual \
					-Wstrict-prototypes -Wmissing-prototypes -Wwrite-strings -Werror
DEBUG_FLAGS=-g -O0
RELEASE_FLAGS=-O3 -DNDEBUG
REQS=convert_dsk.cpp lut.hpp debug.h nanotime.h io.o sort.hpp kmer.hpp utility.hpp
COMPILE=$(CPP) $(CPP_FLAGS)

default: all

lut.hpp: make_lut.py
		python make_lut.py > lut.hpp

io.o: io.hpp io.cpp
		$(COMPILE) $(RELEASE_FLAGS) -c io.cpp

all: $(REQS)
		$(COMPILE) $(RELEASE_FLAGS) -o convert_dsk convert_dsk.cpp io.o

debug: $(REQS)
		$(COMPILE) $(DEBUG_FLAGS) -o convert_dsk convert_dsk.cpp io.o

clean:
		rm -rf convert_dsk *.o *.dSYM lut.hpp
