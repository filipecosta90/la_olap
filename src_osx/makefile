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

grammar.tab.c: grammar.y 
	bison -d grammar.y

dsl: grammar.tab.c lex.yy.c olap.cpp olap.hpp
	$(CC) -o $(EXECNAME) $(CFLAGS) grammar.tab.c lex.yy.c 

.PHONY: clean
clean:
	rm  *.o dsl grammar.tab.c y.output grammar.tab.h grammar.output lex.yy.c dsl  
	
