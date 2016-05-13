
CC = icc
LIBLDIR=
#include
INCLUDE=  
#flags
CFLAGS= -O3 -std=c99 -I${MKLROOT}/include -Iinclude -Ilib/glib-2.0 -Ilib/glib-2.0/include -lglib-2.0 -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel -lmkl_core -lmkl_sequential -lpthread -lm -ldl

querie1: olap.o querie1.o
	$(CC) $(INCLUDE) $(CFLAGS) olap.c querie1.c -o querie1 

all: querie1

.PHONY: clean
clean:
	rm  *.o && rm querie1 && rm -r querie1.dSYM
