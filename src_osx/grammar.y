%{
    
/*************************************************************************
Compiler for the LA language
***************************************************************************/

/*=========================================================================
 C Libraries, Symbol Table, Code Generator & other C code
 =========================================================================*/


#include <stdio.h>
#include <glib.h>
#include <stdlib.h>
#include <string.h>
#include "olap_search.h"
#include "timer.h"
//sleep
#include <unistd.h>

int yylex();
  void yyerror(const char *s);

  double start, stop, elapsed;

%}

%union {
  int ival;
  float fval;
  char *sval;
}

%token BGN END 
%token <sval> IDENTIFIER 
%token <ival> INTEGER
%token HADAMARD KRAO KRON TR
%token VECTOR MATRIX BITMAP
%token BANG TBL_READ MX_FILTER_AND TBL_WRITE CONDITION KEY_CONDITION START STOP

%%

initial_expression : BGN body END
                   ;

body : body elem
     | elem 
     ;

elem : element_declaration ';'
     | time query time
     | atribuition_function ';'
     | function ';'
     ;

query : query atribuition_function ';'
      | query atribuition_expression ';'
      | atribuition_function ';'
      | atribuition_expression ';'
      ;

time : START {
        GET_TIME(start);
      }
     | STOP {
        GET_TIME(stop);
        elapsed = stop - start;
        printf("Tempo: %lf \n",elapsed);
      }
     ;

element_declaration : VECTOR IDENTIFIER
                    | MATRIX IDENTIFIER
                    | BITMAP IDENTIFIER
                    ;

atribuition_function : IDENTIFIER '=' function 
                     ;

atribuition_expression : IDENTIFIER '=' expression
                       ;

expression : IDENTIFIER '*' IDENTIFIER 
           | IDENTIFIER HADAMARD IDENTIFIER 
           | IDENTIFIER KRON IDENTIFIER 
           | IDENTIFIER KRAO IDENTIFIER 
           | IDENTIFIER TR 
           | '(' expression ')'
           ;

function : TBL_READ '(' IDENTIFIER ',' INTEGER ')' {
          char identifier[strlen($3)];
          strcpy(identifier,$3);
          printf("tbl_read(%s,%d)\n",identifier,$5);
}
         | TBL_WRITE '(' IDENTIFIER ',' IDENTIFIER ')'
         | MX_FILTER_AND '(' IDENTIFIER ',' CONDITION ',' KEY_CONDITION ','CONDITION ',' KEY_CONDITION')'
         | BANG '(' INTEGER ')' 
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
