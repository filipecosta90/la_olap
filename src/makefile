CC = icc
LIBLDIR= 

#include
INCLUDE= -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include 

#flags
LIB= -lglib-2.0 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -vec-report7

CFLAGS= $(INCLUDE)  -qopenmp -O3 -g -ggdb -Wall -Wextra  -D_BSD_SOURCE -std=c99 
DEBUGFLAGS= -fopenmp -I${MKLROOT}/include -g -ggdb -O3 -Wall -Wextra -D_BSD_SOURCE -std=c99 -xAVX 

v0: olap_search.o querie1_search.o
	$(CC) $(CFLAGS) -DD_DEBUGGING olap_search.c querie1_search.c -o querie1_search $(LIB) 

v1: olap_search.o querie1_search_v1.o
	$(CC) $(INCLUDE) $(CFLAGS) olap_search.c querie1_search_v1.c -o querie1_search_v1 $(LIB) 

v0debug: olap_search.o querie1_search.o
	$(CC) $(DEBUGFLAGS) -DD_DEBUGGING -o querie1_search olap_search.c querie1_search.c $(LIB) 

all: v0 v1 v0debug

.PHONY: clean
clean:
	rm  *.o  querie1_search querie1_search_v1 querie1_search_debugging
