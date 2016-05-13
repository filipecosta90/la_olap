
CC = icc
LIBLDIR=
#include
INCLUDE=  
#flags
LIBS = -lm -lmkl_intel -lmkl_core -lmkl_sequential -lpthread  -ldl
CFLAGS= -O3 -std=c99 -Ilib/glib-2.0/include -lglib-2.0 -I${MKLROOT}/include -Iinclude -Ilib/glib-2.0  -L${MKLROOT}/lib,-rpath,${MKLROOT}/lib 

querie1: olap.o querie1.o
	$(CC) $(INCLUDE) $(CFLAGS) olap.c querie1.c $(LIBS) -o querie1  

all: querie1

.PHONY: clean
clean:
	rm  *.o && rm querie1 && rm -r querie1.dSYM
