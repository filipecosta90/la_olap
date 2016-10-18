/* The C++ deterministic parser is selected using the skeleton directive */
%skeleton "lalr1.cc"
%require  "3.0"
%debug
%defines
%define api.namespace {OLAP}
%define parser_class_name {OLAP_Parser}

%code requires
{
namespace OLAP {
class OLAP_Driver;
class OLAP_Scanner;
}

// The following definitions is missing when %locations isn't used
# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

}

%parse-param { OLAP_Scanner  &scanner  }
%parse-param { OLAP_Driver  &driver  }

%code{
#include <iostream>
#include <cstdlib>
#include <fstream>

/* include for all driver functions */
#include "olap_driver.hh"

#undef yylex
#define yylex scanner.yylex
}

%define api.value.type variant
%define parse.assert

/* Entry point of grammar */
%start initial_expression

/* Tokens */
%token BGN END 
%token CREATE CUBE 
%token LOAD DROP COLUMN INFILE AS INTO
%token <std::string> IDENTIFIER
%token <int> INTEGER
%token HADAMARD KRAO KRON TR
%token VECTOR MATRIX BITMAP
%token BANG TBL_READ MX_FILTER_AND TBL_WRITE CONDITION KEY_CONDITION START STOP

%locations

%%

initial_expression : BGN body END
                   ;

body : body elem
     | elem
     |
     ;

elem : Create_declaration ';' 
     | Load_declaration ';'
     | Inline_declaration_type ';'
     | time query time
     ;

Create_declaration : CREATE CUBE IDENTIFIER {
                   std::cout << "created cube " << $3 << std::endl;
}
                   ;

Load_declaration : LOAD MATRIX COLUMN INTEGER INFILE IDENTIFIER AS IDENTIFIER INTO IDENTIFIER {
                 driver.load_matrix_csc( $6, $4);
} 
                 | LOAD BITMAP COLUMN INTEGER INFILE IDENTIFIER AS IDENTIFIER INTO IDENTIFIER {
 }
                 ;

Inline_declaration_type : VECTOR IDENTIFIER
                        {
}
| MATRIX IDENTIFIER
{
}
                    | BITMAP IDENTIFIER
                    ;

query : query atribuition_expression ';'
      | atribuition_expression ';'
      ;

time : START {
     }
     | STOP {

}
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

%%

void OLAP::OLAP_Parser::error( const location_type &l, const std::string &err_message )
{
std::cerr << "Error: " << err_message << " at " << l << "\n";
}

