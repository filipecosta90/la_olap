
%union{ char* string; }

%token DOT_P TRANSP  KRAO KRON HADAMARD
%left DOT_P TRANSP KRAO KRON HADAMARD
%token <string> operand_id

%start expression

%%


expression : '(' expression ')'
           | expression DOT_P expression
           | expression TRANSP expression
           | expression KRAO expression
           | expression KRON expression
           | expression HADAMARD expression
           | operand_id
           { printf("%s\n",$1); }
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
