
CC = icc
LIBLDIR=
#include
INCLUDE=  
#flags
LIBS = -lm -lmkl_intel -lmkl_core -lmkl_sequential -lpthread -ldl  -vec-report7
CFLAGS = -O3 -std=c99 -I${MKLROOT}/include -Iinclude 

search: olap_search.o querie1_search.o
	$(CC) $(INCLUDE) $(CFLAGS) olap_search.c querie1_search.c $(LIBS) -o querie1_search


all: search

.PHONY: clean
clean:
	rm  *.o && rm querie1_search && rm -r querie1_search.dSYM
