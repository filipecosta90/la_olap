%token IDENTIFIER CONSTANT STRING_LITERAL 
%token VECTOR MATRIX BITMAP
%token KRAO KRON

%start translation_unit

%%

primary_expression
          : IDENTIFIER
          | CONSTANT
          | STRING_LITERAL
          | '(' expression ')'
          ;

unary_operator
          : 'TR'
          ;

postfix_expression
          : primary_expression
          ;

unary_expression 
          : postfix_expression
          ;

cast_expression
          : unary_expression
          ;

multiplicative_expression
          : cast_expression
          | multiplicative_expression '.' cast_expression
          | multiplicative_expression '><' cast_expression
          | multiplicative_expression KRON cast_expression
          | multiplicative_expression KRAO cast_expression
          ;

additive_expression 
          : multiplicative_expression
          | additive_expression '+' multiplicative_expression
          | additive_expression '-' multiplicative_expression
          ;

assignment_expression
          : unary_expression assignment_operator assignment_expression
          ;

assignment_operator
          : '='
          ;

expression
          : assignment_expression
          | expression ',' assignment_expression
          ;


declaration
          : declaration_specifiers ';'
          | declaration_specifiers init_declarator_list ';'
          ;

declaration_specifiers
          : type_specifier
          | type_specifier declaration_specifiers
          ;

init_declarator_list
          : init_declarator
          | init_declarator_list ',' init_declarator
          ;

init_declarator
          : declarator
          | declarator '=' initializer
          ;

initializer
          : assignment_expression
          ;


type_specifier
          : VECTOR
          | MATRIX
          | BITMAP
          ;

declarator
          : direct_declarator
          ;

direct_declarator
          : IDENTIFIER
          | '(' declarator ')'
          | direct_declarator '(' ')'
          ;

translation_unit
          : external_declaration
          | translation_unit external_declaration
          ;

external_declaration
          : declaration
          ;

%%
#include <stdio.h>

extern char yytext[];
extern int column;

yyerror(s)
char *s;
{
  fflush(stdout);
  printf("\n%*s\n%*s\n", column, "^", column, s);
}
