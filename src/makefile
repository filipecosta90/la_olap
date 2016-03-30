
CC = icc
LIBLDIR=
#include
INCLUDE= 
#flags
CFLAGS= -std=c99 -m64 -Iinclude -Ilib/glib-2.0 -Ilib/glib-2.0/include -lglib-2.0 -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

test: quark_test.o 
	$(CC) $(INCLUDE)  $(CFLAGS) -o quark_test quark_test.c

test.mk: /opt/intel/mkl/bin/mklvars.sh
	    ./$< > $@ intel64

include test.mk

all: test 

.PHONY: clean
clean:
	rm  *.o && rm quark_test mkl
