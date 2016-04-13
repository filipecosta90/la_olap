
CC = icc
LIBLDIR=
#include
INCLUDE= 
#flags
CFLAGS=  -DMKL_ILP64 -m64 -std=c99 -Iinclude -Ilib/glib-2.0 -Ilib/glib-2.0/include -lglib-2.0 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

test: olap.o 
	$(CC) olap.c $(INCLUDE)  $(CFLAGS) -o olap 

test.mk: /opt/intel/mkl/bin/mklvars.sh
	    ./$< > $@ intel64

include test.mk

all: test 

.PHONY: clean
clean:
	rm  *.o && rm olap
