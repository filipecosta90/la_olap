%{
  #include <stdio.h>
  #include <stdlib.h>
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

elem : matrix_declaration ';'
     | atribuition ';'
     | function ';'
     ;

matrix_declaration : type dim idList {printf("Matrix declaration\n");}
                   ;

dim :
    | '(' INTEGER ',' INTEGER ')'
    ;

idList : IDENTIFIER 
       | idList ',' IDENTIFIER
       ;

type : VECTOR
     | MATRIX
     | BITMAP
     ;

atribuition : IDENTIFIER '=' function {printf("Function\n");}
            | IDENTIFIER '=' expression {printf("Expression\n");}
            ;

expression : IDENTIFIER '*' IDENTIFIER {printf("DOT\n");}
           | IDENTIFIER HADAMARD IDENTIFIER {printf("HADAMARD\n");}
           | IDENTIFIER KRON IDENTIFIER {printf("KRON\n");}
           | IDENTIFIER KRAO IDENTIFIER {printf("KRAO\n");}
           | IDENTIFIER TR {printf("TR\n");}
           | '(' expression ')'
           ;

function : IDENTIFIER '(' IDENTIFIER ',' INTEGER ')'
         | IDENTIFIER '(' INTEGER ')'
         | IDENTIFIER '(' IDENTIFIER ')'
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
/*  int s;
  while(s = yylex()){
    printf("%d\n",s);
  }*/
  return (yyparse());
}
