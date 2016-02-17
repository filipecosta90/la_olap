%{
#include <stdio.h>
#include <strings.h>
/* Declaracoes C diversas */

char* expression1, expression2 ;

%}

/* defined matrix as a unidimensional array of integers */

%union{ int* mx; char* string; }

%token <string> id
%type <string> Expression
%%


Expression : '(' Expression ')' 
           { $$ = $2; }
           | '(' Expression ')' 'mult' '(' Expression ')'
           {  expression1 = $2; expression2 = $6; printf(" multiplica\n");}
           | id
           ;

%%
int yyerror(char *s)
{
  fprintf(stderr, "ERRO: %s \n", s);
}

int main()
{
  yyparse();
  return(0);
}
