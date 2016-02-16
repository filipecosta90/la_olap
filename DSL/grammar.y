%{
#include <stdio.h>
#include <strings.h>
/* Declaracoes C diversas */
%}

%token
%type

%%

expression : operation '(' id , id ')'

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
