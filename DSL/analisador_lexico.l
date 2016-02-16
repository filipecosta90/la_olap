%{
#include "y.tab.c" 
/* Declaracoes C diversas */
#define PALAVRA 1
#define PAREN_ESQ 2
#define PONTO 3
#define VIR 4
#define PAREN_DIR 5
%}

%%

[A-Za-z]+ 				{return PALAVRA;}
"("		  				{return PAREN_ESQ;}
"."				 		{return PONTO;}
","						{return VIR;}
")"		  				{return PAREN_DIR;}
.|\n 					{;}

%%
int yywrap(){
    return(1);
}