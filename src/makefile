
CC = icc
LIBLDIR=
#include
INCLUDE=  
#flags
CFLAGS= -g -debug  -I${MKLROOT}/include -Iinclude -Ilib/glib-2.0 -Ilib/glib-2.0/include -lglib-2.0 -std=c99   -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel -lmkl_core -lmkl_sequential -lpthread -lm -ldl


test: olap.o 
	$(CC) $(INCLUDE) $(CFLAGS) olap.c -o olap 

all: test 

.PHONY: clean
clean:
	rm  *.o && rm olap
