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
     IDENTIFIER = 260,
     INTEGER = 261,
     HADAMARD = 262,
     KRAO = 263,
     KRON = 264,
     TR = 265,
     VECTOR = 266,
     MATRIX = 267,
     BITMAP = 268,
     BANG = 269,
     TBL_READ = 270,
     MX_FILTER_AND = 271,
     TBL_WRITE = 272,
     CONDITION = 273,
     KEY_CONDITION = 274,
     START = 275,
     STOP = 276
   };
#endif
/* Tokens.  */
#define BGN 258
#define END 259
#define IDENTIFIER 260
#define INTEGER 261
#define HADAMARD 262
#define KRAO 263
#define KRON 264
#define TR 265
#define VECTOR 266
#define MATRIX 267
#define BITMAP 268
#define BANG 269
#define TBL_READ 270
#define MX_FILTER_AND 271
#define TBL_WRITE 272
#define CONDITION 273
#define KEY_CONDITION 274
#define START 275
#define STOP 276




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 23 "grammar.y"
{
  int ival;
  float fval;
  char *sval;
}
/* Line 1529 of yacc.c.  */
#line 97 "y.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;
