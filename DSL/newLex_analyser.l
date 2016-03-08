D     [0-9]
L     [a-zA-Z_]
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

","     { return(','); }
"="     { return('='); }
"("     { return('('); }
")"     { return(')'); }
"*"     { return('.'); }
"'"     {return (TR);}
"kron"     { return (KRON); }
"><"     { return (HADAMARD); }
"krao"     { return (KRAO); }

[0-9]+  {return (INTEGER);}
L?\"(\\.|[^\\"])*\" { return (IDENTIFIER); }
[ \t\v\n\f]   { }
.     { /* ignore bad characters */ }

%%
