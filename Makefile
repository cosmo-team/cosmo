CC=gcc
C_FLAGS=-m64 -std=c99 -pedantic -W -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual \
				-Wstrict-prototypes -Wmissing-prototypes -Wwrite-strings -Werror
DEBUG_FLAGS=-g
RELEASE_FLAGS= #-O3 -DNDEBUG
COMPILE=$(CC) $(C_FLAGS) $(DEBUG_FLAGS) $(RELEASE_FLAGS)
#COMPILE=$(CC) -o convert_dsk convert_dsk.c io.o transform.o lut.o $(C_FLAGS)

default: all

lut.c: make_lut.py
		python make_lut.py > lut.c

lut.o: lut.h lut.c
		$(COMPILE) -c lut.c

io.o: io.h io.c
		$(COMPILE) -c io.c

sort.o: sort.c sort.h common.h
		$(COMPILE) -c sort.c

transform.o: transform.h transform.c lut.h
		$(COMPILE) -c transform.c

all: convert_dsk.c lut.h debug.h nanotime.h io.o transform.o lut.o sort.o
		$(COMPILE) -o convert_dsk convert_dsk.c io.o transform.o lut.o sort.o

clean:
		rm -rf convert_dsk *.o *.dSYM lut.c
