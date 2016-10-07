/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     BGN = 258,
     END = 259,
     CREATE = 260,
     LOAD = 261,
     DROP = 262,
     IDENTIFIER = 263,
     INTEGER = 264,
     HADAMARD = 265,
     KRAO = 266,
     KRON = 267,
     TR = 268,
     VECTOR = 269,
     MATRIX = 270,
     BITMAP = 271,
     BANG = 272,
     TBL_READ = 273,
     MX_FILTER_AND = 274,
     TBL_WRITE = 275,
     CONDITION = 276,
     KEY_CONDITION = 277,
     START = 278,
     STOP = 279
   };
#endif
/* Tokens.  */
#define BGN 258
#define END 259
#define CREATE 260
#define LOAD 261
#define DROP 262
#define IDENTIFIER 263
#define INTEGER 264
#define HADAMARD 265
#define KRAO 266
#define KRON 267
#define TR 268
#define VECTOR 269
#define MATRIX 270
#define BITMAP 271
#define BANG 272
#define TBL_READ 273
#define MX_FILTER_AND 274
#define TBL_WRITE 275
#define CONDITION 276
#define KEY_CONDITION 277
#define START 278
#define STOP 279




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 28 "grammar.y"
{
  int ival;
  float fval;
  char *sval;
}
/* Line 1529 of yacc.c.  */
#line 103 "y.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

