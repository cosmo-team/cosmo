CC=gcc
C_FLAGS=-m64 -std=c99 -pedantic -W -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual \
				-Wstrict-prototypes -Wmissing-prototypes -Wwrite-strings -Werror
DEBUG_FLAGS=-g
RELEASE_FLAGS=-O3 -DNDEBUG
COMPILE=$(CC) -o convert_dsk convert_dsk.c io.o transform.o lut.o $(C_FLAGS)

default: debug

lut.c: make_lut.py
		python make_lut.py > lut.c

lut.o: lut.h lut.c
		$(CC) $(C_FLAGS) $(DEBUG_FLAGS) -c lut.c

io.o: io.h io.c
		$(CC) $(C_FLAGS) $(DEBUG_FLAGS) -c io.c

transform.o: transform.h transform.c lut.h
		$(CC) $(C_FLAGS) $(DEBUG_FLAGS) -c transform.c

debug: convert_dsk.c lut.h debug.h nanotime.h io.o transform.o lut.o
		$(COMPILE) $(DEBUG_FLAGS)

release: convert_dsk.c lut.h debug.h nanotime.h io.o transform.o lut.o
		$(COMPILE) $(RELEASE_FLAGS)

clean:
		rm -rf convert_dsk *.o *.dSYM lut.c
