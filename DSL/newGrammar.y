%{
  #include <stdio.h>
  #include <stdlib.h>
  #include"lex.yy.c"
  int yylex();
  void yyerror(const char *s);
%}

%token BGN END STRING_LITERAL CHARACTER INTEGER
%token HADAMARD KRAO KRON TR
%token VECTOR MATRIX BITMAP

%%

initial_expression : BGN body END
                   ;

body : matrix_declaration
     | body matrix_declaration
     | atribuition
     ;

matrix_declaration : type identifier
                   | type '(' INTEGER ',' INTEGER ')'
                   ;

type : VECTOR
     | MATRIX
     | BITMAP
     ;

identifier : CHARACTER
           ;

atribuition : identifier '=' function
            | identifier '=' expression
            ;

expression : identifier '.' identifier
           | identifier HADAMARD identifier
           | identifier KRON identifier
           | identifier KRAO identifier
           | identifier TR
           | '(' expression ')'
           ;

function : STRING_LITERAL '(' filename ')'
         ;

filename : STRING_LITERAL
         ;

%%

#include<ctype.h>

void yyerror (const char *s) {
  fprintf (stderr, "%s\n", s);
}

int yywrap(){
  return 1;
}

int main(int argc, char *argv[]){
  yylex();
  return 0;
}
