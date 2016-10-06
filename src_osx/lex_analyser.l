D     [0-9]
H     [a-fA-F0-9]
%{
  #include "y.tab.h"
  #include <string.h>
%}

%%
"START"   {return (START);}
"STOP"    {return (STOP);}
"BEGIN"   {return (BGN);}
"END"     {return (END);}
"vector"  {return (VECTOR);}
"matrix"  {return (MATRIX);}
"bitmap"  {return (BITMAP);}
"bang"    {return (BANG);}
"tbl_read" {return (TBL_READ);}
"tbl_write" {return (TBL_WRITE);}
"mx_filter_and" {return (MX_FILTER_AND);}

";"     { return(';'); }
","     { return(','); }
"="     { return('='); }
"("     { return('('); }
")"     { return(')'); }
"*"     { return('*'); }
"'"     {return (TR);}
"kron"     { return (KRON); }
"><"     { return (HADAMARD); }
"krao"     { return (KRAO); }

[0-9]+  {
          yylval.ival = atoi(yytext);
          return (INTEGER);
        }
[a-zA-Z_'.0-9/]* { 
          char res[strlen(yytext)+1];
          strcpy(res, yytext);
          yylval.sval = res;
          return (IDENTIFIER);
        }
[<=>]+ { return (CONDITION); }
[\t\r]*\'[^']+\' {return(KEY_CONDITION);}
[ \t\v\n\f]   { }
.     { /* ignore bad characters */ }

%%
