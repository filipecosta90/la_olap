D     [0-9]
L     [a-zA-Z_]
H     [a-fA-F0-9]

%{
#include <stdio.h>
#include "y.tab.h"
%}

%%

"vector"     { return("VECTOR"); }
"matrix"     { return("MATRIX"); }
"bitmap"     { return("BITMAP"); }

L?\"(\\.|[^\\"])*\" { count(); return(STRING_LITERAL); }
{L}({L}|{D})*   { count(); return(check_type()); }


";"     { return(';'); }
"{"     { return('{'); }
"}"     { return('}'); }
","     { return(','); }
"="     { return('='); }
"("     { return('('); }
")"     { return(')'); }
"."     { return('.'); }
"-"     { return('-'); }
"+"     { return('+'); }
"*"     { return('*'); }
"/"     { return('/'); }
"kron"     { return("KRON"); }
"><"     { return("><"); }
"krao"     { return("KRAO"); }

[ \t\v\n\f]   { }
.     { /* ignore bad characters */ }

%%

int yywrap(){
    return(1);
}

int check_type()
{
/*
* pseudo code --- this is what it should check
*
* if (yytext == type_name)
*   return(TYPE_NAME);
*
* return(IDENTIFIER);
*/

/*
* it actually will only return IDENTIFIER
*/

  return(IDENTIFIER);
  }


