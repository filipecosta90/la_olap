%{
#include <stdio.h>
#include <strings.h>
/* Declaracoes C diversas */
int yylex();
int yyerror(char *s);
%}

%union{ int* mx; char* string; }

%token DOT_P TRANSP KRAO KNECKER HADAMAR
%left DOT_P TRANSP KRAO KNECKER HADAMAR

%token <string> id

%%


expression : '(' expression ')'
           | expression DOT_P expression
           | TRANSP expression
           | expression KRAO expression
           | expression KNECKER expression
           | expression HADAMAR expression
           |id
           { printf("%s\n",$1);}
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
