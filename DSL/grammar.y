%{
#include <stdio.h>
#include <strings.h>
/* Declaracoes C diversas */
int yylex();
int yyerror(char *s);
%}

/* defined matrix as a unidimensional array of integers */

%union{ int* mx; char* string; }

%token DOT_P TRANSP KRAO KNECKER HADAMAR
%left DOT_P TRANSP KRAO KNECKER HADAMAR

%token <string> id

%%


expression : '(' expression ')' 
           | expression operation expression 
           |id
           { printf("%s\n",$1);}
           ;

operation : DOT_P
         { printf ("invoca multiplicacao mx\n;");}
         | TRANSP
         { printf ("invoca transposicao mx;\n");}
         | KRAO
         { printf ("invoca khatri-rao mx;\n");}
         | KNECKER
         { printf ("invoca khronecker mx;\n");}
         | HADAMAR
         { printf ("invoca hadamar mx;\n");}
         ;

%%
#include "lex.yy.c"
int yyerror(char *s)
{
  fprintf(stderr, "ERRO: %s \n", s);
}

int main()
{
  yyparse();
  return(0);
}
