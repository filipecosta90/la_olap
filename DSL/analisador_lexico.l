%{
%}

id [a-zA-Z]+

%%
[.;(),]         {return(yytext[0]);}
(DOT_P)         {return(DOT_P);}
(TRANSP)        {return(TRANSP);}
(KRAO)          {return(KRAO);}
(KNECKER)       {return(KNECKER);}
(HADAMAR)       {return(HADAMAR);}
{id}            {yylval.string = strdup(yytext); return (id);}
.|\n            {;}

%%
int yywrap(){
    return(1);
}
