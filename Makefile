CC=gcc
C_FLAGS=-m64 -std=c99 -pedantic -W -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual \
				-Wstrict-prototypes -Wmissing-prototypes -Wwrite-strings -Werror
DEBUG_FLAGS=-g
RELEASE_FLAGS=-O3 -DNDEBUG

COMPILE=$(CC) -o convert_dsk convert_dsk.c $(C_FLAGS)

default: debug

lut.h: make_lut.py
		python make_lut.py > lut.h

debug: convert_dsk.c lut.h debug.h
		$(COMPILE) $(DEBUG_FLAGS)

release: convert_dsk.c lut.h debug.h
		$(COMPILE) $(RELEASE_FLAGS)

clean:
		rm -rf convert_dsk *.o *.dSYM lut.h
