
CC = icc
LIBLDIR= -I${MKLROOT}/include -Iinclude -I/home/a57816/lib -I/home/a57816/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include  -lglib-2.0 


#include
INCLUDE= -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include  -lglib-2.0  -I/home/a57816/lib -I/home/a57816/include
 

#flags
LIBS = -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include  -lglib-2.0  -lglib-2.0 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -vec-report7

CFLAGS = -O3 -openmp -std=c99 -farray-notation -xAVX -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include  -lglib-2.0 

search: olap_search.o querie1_search.o
	$(CC) $(INCLUDE) olap_search.c querie1_search.c -o querie1_search $(LIBS) $(CFLAGS) 

all: search

.PHONY: clean
clean:
	rm  *.o && rm querie1_search && rm -r querie1_search.dSYM
