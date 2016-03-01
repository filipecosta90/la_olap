D     [0-9]
L     [a-zA-Z_]
H     [a-fA-F0-9]

%%

"vector"     { return("VECTOR"); }
"matrix"     { return("MATRIX"); }
"bitmap"     { return("BITMAP"); }

L?\"(\\.|[^\\"])*\" { return(STRING_LITERAL); }
{L}({L}|{D})*   {  return(check_type()); }


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
"kron"     { return KRON ; }
"><"     { return("><"); }
"krao"     { return KRAO ; }

[ \t\v\n\f]   { }
.     { /* ignore bad characters */ }

%%

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


