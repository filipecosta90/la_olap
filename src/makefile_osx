
CC = icc
LIBLDIR=
#include
INCLUDE=  
#flags
LIBS = -lm -lmkl_intel -lmkl_core -lmkl_sequential -lpthread -ldl  -vec-report7
CFLAGS = -O3 -std=c99 -Ilib/glib-2.0/include -lglib-2.0 -I${MKLROOT}/include -Iinclude -Ilib/glib-2.0  

osx: olap_osx.o querie1_osx.o
	$(CC) $(INCLUDE) $(CFLAGS) olap_osx.c querie1_osx.c $(LIBS) -o querie1_osx

data: olap_osx.o querie1_data.o
	$(CC) $(INCLUDE) $(CFLAGS) olap_osx.c querie1_data.c $(LIBS) -o querie1_data

all: osx data

.PHONY: clean
clean:
	rm  *.o && rm querie1_osx && rm querie1_data && rm -r querie1.dSYM
