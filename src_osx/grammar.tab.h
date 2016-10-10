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
     CUBE = 261,
     LOAD = 262,
     DROP = 263,
     COLUMN = 264,
     INFILE = 265,
     AS = 266,
     INTO = 267,
     IDENTIFIER = 268,
     INTEGER = 269,
     HADAMARD = 270,
     KRAO = 271,
     KRON = 272,
     TR = 273,
     VECTOR = 274,
     MATRIX = 275,
     BITMAP = 276,
     BANG = 277,
     TBL_READ = 278,
     MX_FILTER_AND = 279,
     TBL_WRITE = 280,
     CONDITION = 281,
     KEY_CONDITION = 282,
     START = 283,
     STOP = 284
   };
#endif
/* Tokens.  */
#define BGN 258
#define END 259
#define CREATE 260
#define CUBE 261
#define LOAD 262
#define DROP 263
#define COLUMN 264
#define INFILE 265
#define AS 266
#define INTO 267
#define IDENTIFIER 268
#define INTEGER 269
#define HADAMARD 270
#define KRAO 271
#define KRON 272
#define TR 273
#define VECTOR 274
#define MATRIX 275
#define BITMAP 276
#define BANG 277
#define TBL_READ 278
#define MX_FILTER_AND 279
#define TBL_WRITE 280
#define CONDITION 281
#define KEY_CONDITION 282
#define START 283
#define STOP 284




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 28 "grammar.y"
{
  int ival;
  float fval;
  char *sval;
}
/* Line 1529 of yacc.c.  */
#line 113 "grammar.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

