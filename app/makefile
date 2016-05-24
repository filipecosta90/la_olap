
CC = icc

BINDIR = 
SRCDIR = src/

LIBLDIR=
#include
INCLUDE=  
#flags
LIBS = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -vec-report7
CFLAGS = -O3 -std=c99 -farray-notation -axAVX -I${MKLROOT}/include -Iinclude 

search: $(SRCDIR)olap_search.o $(SRCDIR)querie1_search.o
	$(CC) $(INCLUDE) $(CFLAGS) $(SRCDIR)olap_search.c $(SRCDIR)querie1_search.c $(LIBS) -o $(BINDIR)querie1_search


all: search

.PHONY: clean
clean:
	rm  $(SRCDIR)*.o && rm $(BINDIR)querie1_search && rm -r querie1_search.dSYM
