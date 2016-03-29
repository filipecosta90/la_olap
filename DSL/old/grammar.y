%{
#include <stdio.h>
#include <stdlib.h>
extern FILE *fp;
%}


%token IDENTIFIER CONSTANT STRING_LITERAL 
%token VECTOR MATRIX BITMAP
%token KRAO KRON HADAMARD TR

%start translation_unit

%%
primary_expression : IDENTIFIER
                   | CONSTANT
                   | STRING_LITERAL
                   | '(' expression ')'
                   ;

unary_operator : TR
               ;

cast_expression : primary_expression
                ;

multiplicative_expression : cast_expression
                          | multiplicative_expression '.' cast_expression
                          | multiplicative_expression HADAMARD cast_expression
                          | multiplicative_expression KRON cast_expression
                          | multiplicative_expression KRAO cast_expression
                          ;

additive_expression : multiplicative_expression
                    | additive_expression '+' multiplicative_expression
                    | additive_expression '-' multiplicative_expression
                    ;

assignment_expression : primary_expression assignment_operator assignment_expression
                      ;

assignment_operator : '='
                    ;

expression : assignment_expression
           | expression ',' assignment_expression
           ;


declaration : declaration_specifiers ';'
            | declaration_specifiers init_declarator_list ';'
            ;

declaration_specifiers : type_specifier
                       | type_specifier declaration_specifiers
                       ;

init_declarator_list : init_declarator
                     | init_declarator_list ',' init_declarator
                     ;

init_declarator : declarator
                | declarator '=' initializer
                ;

initializer : assignment_expression
            ;


type_specifier : VECTOR
               | MATRIX
               | BITMAP
               ;

declarator : direct_declarator
           ;

direct_declarator : IDENTIFIER
                  | '(' declarator ')'
                  | direct_declarator '(' ')'
                  ;

translation_unit : external_declaration
                 | translation_unit external_declaration
                 ;

external_declaration : declaration
                     ;

%%

#include"lex.yy.c"
#include<ctype.h>
int count=0;
int main(int argc, char *argv[])
{
  yyin = fopen(argv[1], "r");

if(!yyparse())
printf("\nParsing complete\n");
else
printf("\nParsing failed\n");
fclose(yyin);
return 0;
}
yyerror(char *s) {
printf("%d : %s %s\n", yylineno, s, yytext );
}         


