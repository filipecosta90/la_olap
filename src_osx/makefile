CC = g++
LIBLDIR= 

#include
INCLUDE= 

#flags
LIB= 
EXECNAME= dsl

CFLAGS= $(INCLUDE)  -Wall -Wextra -O3 

all: dsl 

lex.yy.c: lex_analyser.l
	flex lex_analyser.l

y.tab.c: grammar.y 
	yacc -d grammar.y

dsl: y.tab.c lex.yy.c olap_search.o
	$(CC) -o $(EXECNAME) $(CFLAGS) y.tab.c lex.yy.c 

.PHONY: clean
clean:
	rm  *.o dsl y.tab.c y.output y.tab.h lex.yy.c dsl  
	
