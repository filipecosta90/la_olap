%{
/* C++ string header, for string ops below */
#include <string>

/* Implementation of yyFlexScanner */
#include "olap_parser.hh"
#include "olap_scanner.hh"
#include "olap_driver.hh"

#undef  YY_DECL
#define YY_DECL int OLAP::OLAP_Scanner::yylex( OLAP::OLAP_Parser::semantic_type * const lval, OLAP::OLAP_Parser::location_type *loc )

/* typedef to make the returns for the tokens shorter */
using token = OLAP::OLAP_Parser::token;

/* define yyterminate as this instead of NULL 
#define yyterminate() return( token::END )
*/
/* msvc2010 requires that we exclude this header file. */
#define YY_NO_UNISTD_H

/* update location on matching */
#define YY_USER_ACTION loc->step(); loc->columns(yyleng);

%}

%option debug
%option nodefault
%option yyclass="OLAP::OLAP_Scanner"
%option noyywrap
%option c++

%%
%{          /** Code executed at the beginning of yylex **/
yylval = lval;
%}

 /* The rules.  */

"START"             { return (token::START); }
"STOP"              { return (token::STOP); }
"BEGIN"             { return (token::BGN); }
"END"               { return (token::END); }
"vector"            { return (token::VECTOR); }
"matrix"            { return (token::MATRIX); }
"bitmap"            { return (token::BITMAP); }
"load"              { return (token::LOAD); }
"create"            { return (token::CREATE); }
"cube"              { return (token::CUBE); }
"column"            { return (token::COLUMN); }
"infile"            { return (token::INFILE); }
"as"                { return (token::AS); }
"into"              { return (token::INTO); }
"bang"              { return (token::BANG); }
"tbl_read"          { return (token::TBL_READ); }
"tbl_write"         { return (token::TBL_WRITE); }
"mx_filter_and"     { return (token::MX_FILTER_AND); }
";"                 { return (';'); }
","                 { return (','); }
"="                 { return ('='); }
"("                 { return ('('); }
")"                 { return (')'); }
"*"                 { return ('*'); }
"'"                 { return (token::TR); }
"kron"              { return (token::KRON); }
"><"                { return (token::HADAMARD); }
"krao"              { return (token::KRAO); }

[0-9]+              {
      int i_auto = std::stoi (yytext,nullptr,0);
      yylval->build<int>(i_auto);
      return (token::INTEGER);
}

[a-zA-Z_'.0-9/]*    { 
  yylval->build< std::string >( yytext );
  return (token::IDENTIFIER);
}

[<=>]+              { return (token::CONDITION); }
[\t\r]*\'[^']+\'    {return(token::KEY_CONDITION);}
[ \t\v\n\f]         { /* ignore spacing */ }
.                   { /* ignore bad characters */ }

%%

