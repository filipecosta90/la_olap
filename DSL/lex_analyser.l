%{
%}

operand_id [a-zA-Z]+

%%
[.;(),]  {return(yytext[0]);}
"dot"   {return(DOT_P);}
"'"  {return(TRANSP);}
"krao"   {return(KRAO);}
"kron"   {return(KRON);}
"><" {return(HADAMARD);}
{operand_id} {yylval.string = strdup(yytext); return (operand_id);}
.|\n            {;}

%%
int yywrap(){
    return(1);
}
