CC=gcc
C_FLAGS=-m64 -std=c99 -pedantic -W -Wall -Wextra -Wshadow -Wpointer-arith -Wcast-qual \
				-Wstrict-prototypes -Wmissing-prototypes -Wwrite-strings -Werror
RELEASE_FLAGS=-O3 -DNDEBUG

COMPILE=$(CC) -o convert_dsk convert_dsk.c $(C_FLAGS)

default: debug

debug: convert_dsk.c debug.h
		$(COMPILE)

release: convert_dsk.c debug.h
		$(COMPILE) $(RELEASE_FLAGS)

clean:
		rm -rf convert_dsk *.o
