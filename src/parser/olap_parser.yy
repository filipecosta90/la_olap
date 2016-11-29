/* The C++ deterministic parser is selected using the skeleton directive */
%skeleton "lalr1.cc"
%require  "3.0"
%debug
%defines
%define api.namespace {OLAP}
%define parser_class_name {OLAP_Parser}

%{

#include <glib.h>
#include "olap_parser.hh"
#include "olap_scanner.hh"
#include "olap_engine.hxx"

#undef yylex
#define yylex scanner.yylex

// The following definitions is missing when %locations isn't used
# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

OLAP::OLAP_Engine *db_engine = NULL;
GHashTable *cubes_table = NULL; 
%}

%code requires
{

#include <iostream>
#include <cstdlib>
#include <fstream>

/* include for all driver functions */
#include "olap_driver.hh"

}

%code provides
{
namespace OLAP {
// Forward declaration of the Driver and Scanner classes
class OLAP_Driver;
class OLAP_Parser;
class OLAP_Engine;
}
}
%parse-param { OLAP_Scanner  &scanner  }
%parse-param { OLAP_Driver  &driver  }

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

initial_expression : BGN 
                        {    
                        db_engine = new OLAP::OLAP_Engine();
                       }
                        body END
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
  db_engine->create_cube($3);
 std::cout << "created cube " << $3 << std::endl;
   }
                   ;

Load_declaration : LOAD MATRIX COLUMN INTEGER INFILE IDENTIFIER AS IDENTIFIER INTO IDENTIFIER {                 
                 OLAP::OLAP_Cube *cube = db_engine->cube_lookup($10);
                 cube->load_matrix_to_csc_from_tbl($6,$4,100);
                 std::cout << "load  into " << $10 << std::endl;
               } 
                 | LOAD BITMAP COLUMN INTEGER INFILE IDENTIFIER AS IDENTIFIER INTO IDENTIFIER {
}
                 ;

Inline_declaration_type : VECTOR IDENTIFIER {
                        }
| MATRIX IDENTIFIER {
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

