
CC = icc
LIBLDIR=
#include
INCLUDE=  
#flags
CFLAGS= -g -debug  -I${MKLROOT}/include -Iinclude -Ilib/glib-2.0 -Ilib/glib-2.0/include -lglib-2.0 -std=c99   -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel -lmkl_core -lmkl_sequential -lpthread -lm -ldl


test: olap.o 
	$(CC) $(INCLUDE) $(CFLAGS) olap.c -o olap 

querie1: querie1.o
	$(CC) $(INCLUDE) $(CFLAGS) querie1.c -o querie1 

all: test querie1

.PHONY: clean
clean:
	rm  *.o && rm olap && rm querie1
