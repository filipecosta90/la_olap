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
%token CREATE LOAD DROP 
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

elem : Inline_declaration_type ';'
     | time query time
     | atribuition_function ';'
     | function ';'
     ;

Inline_declaration_type : VECTOR IDENTIFIER
{
   printf("need to create vector variable\n");
   //CSC
    __declspec(align(MEM_LINE_SIZE)) float* vector_csc_values;
    __declspec(align(MEM_LINE_SIZE)) int* vector_row_ind;
    //COMMON
    int vector_n_nnz;
    int vector_n_rows;
}
| MATRIX IDENTIFIER
{
   printf("need to create matrix variable\n");
    //CSC
    __declspec(align(MEM_LINE_SIZE)) float* matrix_csc_values;
    __declspec(align(MEM_LINE_SIZE)) int* matrix_row_ind;
    __declspec(align(MEM_LINE_SIZE)) int* matrix_col_ptr;
    //COMMON
    int matrix_n_nnz;
    int matrix_n_rows;
    int matrix_n_cols;
}
                    | BITMAP IDENTIFIER
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
