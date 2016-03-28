
CC = gcc-5

#include
INCLUDE= -I/usr/local/include/glib-2.0 -I/usr/local/lib/glib-2.0/include 

#flags
CFLAGS= -m64 -I/usr/local/include/glib-2.0 -I/usr/local/lib/glib-2.0/include -lglib-2.0 -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

test: quark_test.o 
	$(CC) $(INCLUDE) $(CFLAGS) -o quark_test quark_test.c

mkl: mkl-lab-solution.o 
	$(CC) $(INCLUDE) $(CFLAGS) -o mkl mkl-lab-solution.c


all: test zeros

.PHONY: clean
clean:
	rm  *.o && rm quark_test mkl
