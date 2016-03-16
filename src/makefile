
CC = gcc-5

#include
INCLUDE= -I/usr/local/include/glib-2.0 -I/usr/local/lib/glib-2.0/include 

#flags
CFLAGS= -I/usr/local/include/glib-2.0 -I/usr/local/lib/glib-2.0/include -lglib-2.0 

quark_test: quark_test.o 
	$(CC) $(INCLUDE) $(CFLAGS) -o quark_test quark_test.c

all: trapezoid create

.PHONY: clean
clean:
	rm  *.o && rm quark_test 
