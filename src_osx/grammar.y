%{
    
/*************************************************************************
Compiler for the LA language
***************************************************************************/

/*=========================================================================
 C Libraries, Symbol Table, Code Generator & other C code
 =========================================================================*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "olap_search.hpp"
#include "timer.hpp"
//sleep:
#include <unistd.h>

using namespace std;
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
%token CREATE CUBE 
%token LOAD DROP COLUMN INFILE AS INTO
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

elem : Create_declaration ';' 
     | Load_declaration ';'
     | Inline_declaration_type ';'
     | time query time
     | atribuition_function ';'
     | function ';'
     ;

Create_declaration : CREATE CUBE IDENTIFIER {} 
                   ;

Load_declaration : LOAD MATRIX COLUMN INTEGER INFILE IDENTIFIER AS IDENTIFIER INTO IDENTIFIER {} 
                 | LOAD BITMAP COLUMN INTEGER INFILE IDENTIFIER AS IDENTIFIER INTO IDENTIFIER { }
                 ;

Inline_declaration_type : VECTOR IDENTIFIER
{
   //CSC
/*    __declspec(align(MEM_LINE_SIZE)) float* vector_csc_values;
    __declspec(align(MEM_LINE_SIZE)) int* vector_row_ind;
    //COMMON
    int vector_n_nnz;
    int vector_n_rows;
*/
}
| MATRIX IDENTIFIER
{
/*
printf("need to create matrix variable\n");
    //CSC
    __declspec(align(MEM_LINE_SIZE)) float* matrix_csc_values;
    __declspec(align(MEM_LINE_SIZE)) int* matrix_row_ind;
    __declspec(align(MEM_LINE_SIZE)) int* matrix_col_ptr;
    //COMMON
    int matrix_n_nnz;
    int matrix_n_rows;
    int matrix_n_cols;
*/
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


void yyerror (const char *s) {
  fprintf (stderr, "%s\n", s);
}

int yywrap() {
                        return 1;
                        }


int main(int argc, char *argv[]){
/*  int s;
  while(s = yylex()){
    printf("%d\n",s);
  }*/
  return (yyparse());
}
