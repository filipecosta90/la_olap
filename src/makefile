CC = icc
LIBLDIR= 

#include
INCLUDE= -I${MKLROOT}/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include 

#flags
LIB= -lglib-2.0 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm 

CFLAGS=  -openmp  -O3 -Wall -Wextra -D_BSD_SOURCE -std=c99 -xAVX $(INCLUDE) -vec-report7
DEBUGFLAGS= -fopenmp -I${MKLROOT}/include -g -ggdb -O3 -Wall -Wextra -D_BSD_SOURCE -std=c99 -xAVX $(INCLUDE) 

v0: olap_search.o querie1_search.o
	$(CC) $(CFLAGS) olap_search.c querie1_search.c -o querie1_search $(LIB) 

v1: olap_search.o querie1_search_v1.o
	$(CC) $(INCLUDE) $(CFLAGS) olap_search.c querie1_search_v1.c -o querie1_search_v1 $(LIB) 

v0debug: olap_search.o querie1_search.o
	$(CC) $(DEBUGFLAGS) -DD_DEBUGGING -o querie1_search_debugging olap_search.c querie1_search.c $(LIB) 

all: v0 v1 v0debug

.PHONY: clean
clean:
	rm  *.o  querie1_search querie1_search_v1 querie1_search_debugging
