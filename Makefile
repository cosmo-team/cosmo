CC=gcc
C_FLAGS=-m64 -std=c99 -pedantic -W -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual \
				-Wstrict-prototypes -Wmissing-prototypes -Wwrite-strings -Werror
DEBUG_FLAGS=-g
RELEASE_FLAGS=-O3 -DNDEBUG
EXPERIMENT_FLAGS=-DEXPERIMENT
COMPILE=$(CC) -o convert_dsk convert_dsk.c io.o $(C_FLAGS)

default: debug

lut.h: make_lut.py
		python make_lut.py > lut.h

io.o: io.h io.c
		$(CC) $(C_FLAGS) $(DEBUG_FLAGS) -c io.c

debug: convert_dsk.c lut.h debug.h nanotime.h io.o
		$(COMPILE) $(DEBUG_FLAGS)

release: convert_dsk.c lut.h debug.h nanotime.h io.o
		$(COMPILE) $(RELEASE_FLAGS)

clean:
		rm -rf convert_dsk *.o *.dSYM lut.h
