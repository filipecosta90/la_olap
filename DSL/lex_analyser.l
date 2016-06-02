D     [0-9]
H     [a-fA-F0-9]
%{
  #include "y.tab.h"
%}

%%
"BEGIN"   {return (BGN);}
"END"     {return (END);}
"vector"  {return (VECTOR);}
"matrix"  {return (MATRIX);}
"bitmap"  {return (BITMAP);}
"bang"    {return (BANG);}
"tbl_read" {return (TBL_READ);}
"tbl_write" {return (TBL_WRITE);}
"tbl_filter" {return (TBL_FILTER);}

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

[0-9]+  {return (INTEGER);}
[a-zA-Z_'.0-9/]* { return (IDENTIFIER); }
[<=>]+ { return (CONDITION); }
[\t\r]*\'[^']+\' {return(KEY_CONDITION);}
[ \t\v\n\f]   { }
.     { /* ignore bad characters */ }

%%

/*int yywrap(){return 1;}

int main(){
  int s;
  while( s = yylex() ){
  printf("%d\n",s);
  }
}*/
