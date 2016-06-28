CC = icc
LIBLDIR= 

#include
INCLUDE= -I${MKLROOT}/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include 

#flags
LIB= -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -lglib-2.0 

CFLAGS= -qopenmp -I${MKLROOT}/include -g -ggdb -O3 -Wall -Wextra -D_BSD_SOURCE -std=c99 -farray-notation -xAVX $(INCLUDE) 

v0: olap_search.o querie1_search.o
	$(CC) $(CFLAGS) olap_search.c querie1_search.c -o querie1_search $(LIB) 

v1: olap_search.o querie1_search_v1.o
	$(CC) $(INCLUDE) $(CFLAGS) olap_search.c querie1_search_v1.c -o querie1_search_v1 $(LIB) 

v2: olap_search.o querie1_search.o
	$(CC) $(CFLAGS) -DD_DEBUGGING olap_search.c querie1_search.c -o querie1_search_debugging $(LIB) 

all: v0 v1 v2

.PHONY: clean
clean:
	rm  *.o  querie1_search querie1_search_v1 querie1_search_debugging
