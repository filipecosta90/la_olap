CC = icc
LIBLDIR= 

#include
INCLUDE= `pkg-config --cflags --libs glib-2.0` 

#flags
LIB= 
EXECNAME= dsl
SEQCFLAGS= $(INCLUDE) -O3 -Wall -Wextra -D_BSD_SOURCE -std=c99 
CFLAGS= $(INCLUDE)  -fopenmp -Wall -Wextra -D_BSD_SOURCE -std=c99 
DEBUGFLAGS=  $(INCLUDE) -fopenmp -g -ggdb -O3 -Wall -Wextra -D_BSD_SOURCE -std=c99 

all: dsl 

y.tab.c: lex_analyser.l
	flex lex_analyser.l

lex.yy.c: grammar.y y.tab.c
	yacc -d -v grammar.y

dsl: y.tab.c lex.yy.c olap_search.o
	$(CC) -o $(EXECNAME) $(CFLAGS) y.tab.c lex.yy.c

nomkl: olap_search.o query1_no_mkl.o
	$(CC) $(CFLAGS) olap_search.c query1_no_mkl.c -o query1_no_mkl $(LIB) 

.PHONY: clean
clean:
	rm  *.o dsl y.tab.c y.output y.tab.h lex.yy.c algebric
	
