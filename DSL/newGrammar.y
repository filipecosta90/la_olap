%{
  #include <stdio.h>
  #include <stdlib.h>
  #include"lex.yy.c"
  int yylex();
  void yyerror(const char *s);
%}

%token BGN END IDENTIFIER INTEGER
%token HADAMARD KRAO KRON TR
%token VECTOR MATRIX BITMAP

%%

initial_expression : BGN body END
                   ;

body : body elem
     | elem 
     ;

elem : matrix_declaration
     | atribuition
     | function
     ;

matrix_declaration : type dim idList
                   ;
dim  : 
     | '(' INTEGER ',' INTEGER ')'
     ;

idList : IDENTIFIER 
       | idList ',' IDENTIFIER
       ;

type : VECTOR
     | MATRIX
     | BITMAP
     ;

atribuition : IDENTIFIER '=' function
            | IDENTIFIER '=' expression
            ;

expression : IDENTIFIER '*' IDENTIFIER
           | IDENTIFIER HADAMARD IDENTIFIER
           | IDENTIFIER KRON IDENTIFIER
           | IDENTIFIER KRAO IDENTIFIER
           | IDENTIFIER TR
           | '(' expression ')'
           ;

function : IDENTIFIER '(' IDENTIFIER ')'
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
