CXX=g++
CPP_FLAGS=-m64 -std=c++0x -pedantic-errors -W -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual \
					-Wunused -Wstrict-prototypes -Wmissing-prototypes -Wwrite-strings -Werror
DEBUG_FLAGS=-g -O0
# The MMX and SSE flags can be safely disabled
RELEASE_FLAGS=-O3 -DNDEBUG -mmmx -msse -msse2 -msse3 -msse4 -march=native
REQS=kramer.cpp lut.hpp debug.h nanotime.h io.o sort.hpp kmer.hpp dummies.hpp
COMPILE=$(CXX) $(CPP_FLAGS)

default: all

lut.hpp: make_lut.py
		python make_lut.py > lut.hpp

io.o: io.hpp io.cpp debug.h dummies.hpp kmer.hpp
		$(COMPILE) $(RELEASE_FLAGS) -c io.cpp

all: $(REQS)
		$(COMPILE) $(RELEASE_FLAGS) -o kramer kramer.cpp io.o

debug: $(REQS)
		$(COMPILE) $(DEBUG_FLAGS) -o kramer kramer.cpp io.o

clean:
		rm -rf kramer *.o *.dSYM lut.hpp
