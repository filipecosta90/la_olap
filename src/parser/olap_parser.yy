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
                        {    cubes_table = g_hash_table_new_full(
                              g_str_hash, g_str_equal, //< This is an integer hash.
                                    free, //< Call "free" on the key (made with "malloc").
                                          free //< Call "free" on the value (made with "strdup").
                                             );
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
                   GHashTable *cube = g_hash_table_new_full (g_str_hash,
        g_str_equal, g_free, g_free);
          g_hash_table_insert (cubes_table,
                               (void*) $3.c_str() ,
                                                    cube);
std::cout << "created cube " << $3 << std::endl;
   }
                   ;

Load_declaration : LOAD MATRIX COLUMN INTEGER INFILE IDENTIFIER AS IDENTIFIER INTO IDENTIFIER {
                 
                   GHashTable *cube = (GHashTable*) g_hash_table_lookup ( cubes_table,  (void*) $10.c_str() );
                 std::cout << "load  into " << $10 << std::endl;
                 int n_nnz, n_rows, n_cols;
                 float *A_csc_values;
                 int *A_row_ind, *A_col_ptr;
                 driver.load_matrix_csc( $6, $4, 100, &n_nnz, &n_rows, &n_cols, &A_csc_values, &A_row_ind, &A_col_ptr );
                 std::cout << "NNZ " << n_nnz << "|\tN_ROWS " << n_rows << "|\tN_COLS " << n_cols << std::endl;

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

