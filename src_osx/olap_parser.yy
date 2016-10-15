/* C++ parser interface */
%skeleton "lalr1.cc"

/* require bison version */
%require  "3.0"

/* call yylex with a location */
%locations

%debug
%defines
%define api.namespace {OLAP}
%define parser_class_name {OLAP_Parser}
%define api.value.type variant

/* assert correct cleanup of semantic value objects */
%define parse.assert

%parse-param { OLAP_Scanner  &scanner  }
%parse-param { OLAP_Driver  &driver  }

/* inserted near top of header + source file */
%code requires{
namespace OLAP {
class OLAP_Driver;
class OLAP_Scanner;
class OLAP_LA;
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

/* inserted near top of source file */
%code{
#include <iostream>
#include <cstdlib>
#include <fstream>

/* include for all driver functions */
#include "olap_driver.hpp"

#undef yylex
#define yylex scanner.yylex
}

%token BGN END 
%token CREATE CUBE 
%token LOAD DROP COLUMN INFILE AS INTO
%token <std::string> IDENTIFIER
%token <int> INTEGER
%token HADAMARD KRAO KRON TR
%token VECTOR MATRIX BITMAP
%token BANG TBL_READ MX_FILTER_AND TBL_WRITE CONDITION KEY_CONDITION START STOP

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


void
OLAP::OLAP_Parser::error( const location_type &l, const std::string &err_message )
    {
std::cerr << "Error: " << err_message << " at " << l << "\n";
   }

