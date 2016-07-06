CC = icc
LIBLDIR= 

#include
INCLUDE= -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include 

#flags
LIB= -lglib-2.0 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -vec-report7

SEQCFLAGS= $(INCLUDE) -O3 -Wall -Wextra -xAVX -D_BSD_SOURCE -std=c99 
CFLAGS= $(INCLUDE)  -qopenmp -openmp-report7 -O3 -Wall -Wextra -xAVX -D_BSD_SOURCE -std=c99 
DEBUGFLAGS=  $(INCLUDE) -qopenmp -g -ggdb -O3 -Wall -Wextra -xAVX -D_BSD_SOURCE -std=c99 


nomkl: olap_search.o query1_no_mkl.o
	$(CC) $(CFLAGS) olap_search.c query1_no_mkl.c -o query1_no_mkl $(LIB) 

nomklver: olap_search.o query1_no_mkl.o
	$(CC) $(CFLAGS) -DD_VERBOSE olap_search.c query1_no_mkl.c -o query1_no_mkl $(LIB) 

seqnomklver: olap_search.o query1_no_mkl.o
	$(CC) $(SEQCFLAGS) -DD_VERBOSE olap_search.c query1_no_mkl.c -o query1_no_mkl $(LIB) 


nomkldebug: olap_search.o query1_no_mkl.o
	$(CC) $(DEBUGFLAGS) olap_search.c query1_no_mkl.c -o query1_no_mkl $(LIB) 

nomkldebugver: olap_search.o query1_no_mkl.o
	$(CC) $(DEBUGFLAGS) -DD_VERBOSE -DD_DEBUGGING olap_search.c query1_no_mkl.c -o query1_no_mkl $(LIB) 


v0: olap_search.o querie1_search.o
	$(CC) $(CFLAGS) -DD_DEBUGGING olap_search.c querie1_search.c -o querie1_search $(LIB) 

v1: olap_search.o querie1_search_v1.o
	$(CC) $(INCLUDE) $(CFLAGS) olap_search.c querie1_search_v1.c -o querie1_search_v1 $(LIB) 

v0debug: olap_search.o querie1_search.o
	$(CC) $(DEBUGFLAGS) -DD_DEBUGGING -o querie1_search olap_search.c querie1_search.c $(LIB) 

all: v0 v1 v0debug

.PHONY: clean
clean:
	rm  *.o query1_no_mkl querie1_search querie1_search_v1 querie1_search
