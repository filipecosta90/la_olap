/* ---------------------------------------------------------------------------
 **    Filename: olap_scanner.hh
 **
 **     License: This file is part of OLAP PROJECT.
 **
 **              OLAP PROJECT is free software: you can redistribute it
 **              and/or modify it under the terms of the GNU General Public
 **              License as published by the Free Software Foundation,
 **              either version 3 of the License, or (at your option)
 **              any later version.
 **
 **              OLAP is distributed in the hope that it will be useful,
 **              but WITHOUT ANY WARRANTY; without even the implied warranty of
 **              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **              GNU General Public License for more details.
 **
 **              You should have received a copy of the GNU General Public
 **              License along with OLAP.
 **              If not, see <http://www.gnu.org/licenses/>.
 **
 ** Description:
 **
 **     Authors: Filipe Oliveira <a57816@alunos.uminho.pt>
 **
 ** University of Minho, High Performance Computing Dpt. , October 2016
 ** -------------------------------------------------------------------------*/

#ifndef __OLAPSCANNER_HH__
#define __OLAPSCANNER_HH__ 1

#if ! defined(yyFlexLexerOnce)
#include <FlexLexer.h>
#endif

#include "olap_parser.hh"
#include "location.hh"

namespace OLAP
{
  class OLAP_Scanner : public yyFlexLexer {
  public:

      OLAP_Scanner(std::istream *in) : yyFlexLexer(in)
      {
          loc = new OLAP::OLAP_Parser::location_type();
      };
      virtual ~OLAP_Scanner() {
          delete loc;
      };
      
      //get rid of override virtual function warning
      using FlexLexer::yylex;
      
      virtual
      int yylex( OLAP::OLAP_Parser::semantic_type * const lval,
                OLAP::OLAP_Parser::location_type *location );
      // YY_DECL defined in olap_lexer.l
      // Method body created by flex in olap_lexer.yy.cc
      
      
  private:
      /* yyval ptr */
      OLAP::OLAP_Parser::semantic_type *yylval = nullptr;
      /* location ptr */
      OLAP::OLAP_Parser::location_type *loc    = nullptr;
  };
    
} /* end namespace OLAP */

#endif /* END __OLAPSCANNER_HH__ */
