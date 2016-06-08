
CC = icc
LIBLDIR= 


#include
INCLUDE= -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include  -lglib-2.0  
 

#flags
LIB = -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include -lglib-2.0 -lglib-2.0 -lmkl_intel_ilp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl #-vec-report7

CFLAGS = -O3 -Wall -Wextra -qopenmp -D_BSD_SOURCE -std=c99 -farray-notation -xAVX  -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include -lglib-2.0 -lglib-2.0 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread # -vec-report7
search0: olap_search.o querie1_search.o
	$(CC) $(INCLUDE) olap_search.c querie1_search.c -o querie1_search $(CFLAGS) $(LIB) 

search1: olap_search.o querie1_search_v1.o
	$(CC) $(INCLUDE) olap_search.c querie1_search_v1.c -o querie1_search_v1 $(CFLAGS) $(LIB) 


all: search0 search1

.PHONY: clean
clean:
	rm  *.o  querie1_search querie1_search_v1
