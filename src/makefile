
CC = icc
LIBLDIR= 


#include
INCLUDE=  -I${MKLROOT}/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include 

#flags
LIB =  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -lglib-2.0 

CFLAGS = $(INCLUDE) -openmp -g -ggdb -O3 -Wall -Wextra -D_BSD_SOURCE -std=c99 -farray-notation -xAVX  

v0: olap_search.o querie1_search.o
	$(CC) $(INCLUDE) $(CFLAGS) olap_search.c querie1_search.c -o querie1_search $(LIB) 

v1: olap_search.o querie1_search_v1.o
	$(CC) $(INCLUDE) $(CFLAGS) olap_search.c querie1_search_v1.c -o querie1_search_v1 $(LIB) 


all: v0 v1

.PHONY: clean
clean:
	rm  *.o  querie1_search querie1_search_v1
