CC    ?= clang
CXX   ?= clang++

EXE = dsl

CDEBUG = -g -Wall

CXXDEBUG = -g -Wall

CSTD = -std=c99
CXXSTD = -std=c++14

CFLAGS = -Wno-deprecated-register -O3  $(CDEBUG) $(CSTD)
CXXFLAGS = -Wno-deprecated-register -O3  $(CXXDEBUG) $(CXXSTD)


CPPOBJ = main olap_driver
SOBJ =  parser lexer

FILES = $(addsuffix .cpp, $(CPPOBJ))

OBJS  = $(addsuffix .o, $(CPPOBJ))

CLEANLIST =  $(addsuffix .o, $(OBJ)) $(OBJS) \
						 olap_parser.tab.cc olap_parser.tab.hh \
						 location.hh position.hh \
						 stack.hh olap_parser.output parser.o \
						 lexer.o olap_lexer.yy.cc $(EXE)\

.PHONY: all
all: dsl

dsl: $(FILES)
	$(MAKE) $(SOBJ)
	$(MAKE) $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(OBJS) parser.o lexer.o $(LIBS)


parser: olap_parser.yy
	bison -d -v olap_parser.yy
	$(CXX) $(CXXFLAGS) -c -o parser.o olap_parser.tab.cc

lexer: olap_lexer.l
	flex --outfile=olap_lexer.yy.cc  $<
	$(CXX)  $(CXXFLAGS) -c olap_lexer.yy.cc -o lexer.o

.PHONY: clean
clean:
	rm -rf $(CLEANLIST)

