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
#include "olap.hpp"
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

Load_declaration : LOAD MATRIX COLUMN INTEGER INFILE IDENTIFIER AS IDENTIFIER INTO IDENTIFIER {
  col_read_csc ( $6, $4 );
                 
                 } 
                 | LOAD BITMAP COLUMN INTEGER INFILE IDENTIFIER AS IDENTIFIER INTO IDENTIFIER { }
                 ;

Inline_declaration_type : VECTOR IDENTIFIER
{
}
| MATRIX IDENTIFIER
{
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
  return (yyparse());
}
