// A Bison parser, made by GNU Bison 3.0.4.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.


// First part of user declarations.

#line 37 "olap_parser.tab.cc" // lalr1.cc:404

# ifndef YY_NULLPTR
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULLPTR nullptr
#  else
#   define YY_NULLPTR 0
#  endif
# endif

#include "olap_parser.tab.hh"

// User implementation prologue.

#line 51 "olap_parser.tab.cc" // lalr1.cc:412
// Unqualified %code blocks.
#line 43 "olap_parser.yy" // lalr1.cc:413

#include <iostream>
#include <cstdlib>
#include <fstream>

/* include for all driver functions */
#include "olap_driver.hpp"

#undef yylex
#define yylex scanner.yylex

#line 65 "olap_parser.tab.cc" // lalr1.cc:413


#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K].location)
/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

# ifndef YYLLOC_DEFAULT
#  define YYLLOC_DEFAULT(Current, Rhs, N)                               \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).begin  = YYRHSLOC (Rhs, 1).begin;                   \
          (Current).end    = YYRHSLOC (Rhs, N).end;                     \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).begin = (Current).end = YYRHSLOC (Rhs, 0).end;      \
        }                                                               \
    while (/*CONSTCOND*/ false)
# endif


// Suppress unused-variable warnings by "using" E.
#define YYUSE(E) ((void) (E))

// Enable debugging if requested.
#if YYDEBUG

// A pseudo ostream that takes yydebug_ into account.
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Symbol)         \
  do {                                          \
    if (yydebug_)                               \
    {                                           \
      *yycdebug_ << Title << ' ';               \
      yy_print_ (*yycdebug_, Symbol);           \
      *yycdebug_ << std::endl;                  \
    }                                           \
  } while (false)

# define YY_REDUCE_PRINT(Rule)          \
  do {                                  \
    if (yydebug_)                       \
      yy_reduce_print_ (Rule);          \
  } while (false)

# define YY_STACK_PRINT()               \
  do {                                  \
    if (yydebug_)                       \
      yystack_print_ ();                \
  } while (false)

#else // !YYDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YYUSE(Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void>(0)
# define YY_STACK_PRINT()                static_cast<void>(0)

#endif // !YYDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)

#line 12 "olap_parser.yy" // lalr1.cc:479
namespace OLAP {
#line 151 "olap_parser.tab.cc" // lalr1.cc:479

  /// Build a parser object.
  OLAP_Parser::OLAP_Parser (OLAP_Scanner  &scanner_yyarg, OLAP_Driver  &driver_yyarg)
    :
#if YYDEBUG
      yydebug_ (false),
      yycdebug_ (&std::cerr),
#endif
      scanner (scanner_yyarg),
      driver (driver_yyarg)
  {}

  OLAP_Parser::~OLAP_Parser ()
  {}


  /*---------------.
  | Symbol types.  |
  `---------------*/

  inline
  OLAP_Parser::syntax_error::syntax_error (const location_type& l, const std::string& m)
    : std::runtime_error (m)
    , location (l)
  {}

  // basic_symbol.
  template <typename Base>
  inline
  OLAP_Parser::basic_symbol<Base>::basic_symbol ()
    : value ()
  {}

  template <typename Base>
  inline
  OLAP_Parser::basic_symbol<Base>::basic_symbol (const basic_symbol& other)
    : Base (other)
    , value ()
    , location (other.location)
  {
      switch (other.type_get ())
    {
      case 14: // INTEGER
        value.copy< int > (other.value);
        break;

      case 13: // IDENTIFIER
        value.copy< std::string > (other.value);
        break;

      default:
        break;
    }

  }


  template <typename Base>
  inline
  OLAP_Parser::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, const semantic_type& v, const location_type& l)
    : Base (t)
    , value ()
    , location (l)
  {
    (void) v;
      switch (this->type_get ())
    {
      case 14: // INTEGER
        value.copy< int > (v);
        break;

      case 13: // IDENTIFIER
        value.copy< std::string > (v);
        break;

      default:
        break;
    }
}


  // Implementation of basic_symbol constructor for each type.

  template <typename Base>
  OLAP_Parser::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, const location_type& l)
    : Base (t)
    , value ()
    , location (l)
  {}

  template <typename Base>
  OLAP_Parser::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, const int v, const location_type& l)
    : Base (t)
    , value (v)
    , location (l)
  {}

  template <typename Base>
  OLAP_Parser::basic_symbol<Base>::basic_symbol (typename Base::kind_type t, const std::string v, const location_type& l)
    : Base (t)
    , value (v)
    , location (l)
  {}


  template <typename Base>
  inline
  OLAP_Parser::basic_symbol<Base>::~basic_symbol ()
  {
    clear ();
  }

  template <typename Base>
  inline
  void
  OLAP_Parser::basic_symbol<Base>::clear ()
  {
    // User destructor.
    symbol_number_type yytype = this->type_get ();
    basic_symbol<Base>& yysym = *this;
    (void) yysym;
    switch (yytype)
    {
   default:
      break;
    }

    // Type destructor.
    switch (yytype)
    {
      case 14: // INTEGER
        value.template destroy< int > ();
        break;

      case 13: // IDENTIFIER
        value.template destroy< std::string > ();
        break;

      default:
        break;
    }

    Base::clear ();
  }

  template <typename Base>
  inline
  bool
  OLAP_Parser::basic_symbol<Base>::empty () const
  {
    return Base::type_get () == empty_symbol;
  }

  template <typename Base>
  inline
  void
  OLAP_Parser::basic_symbol<Base>::move (basic_symbol& s)
  {
    super_type::move(s);
      switch (this->type_get ())
    {
      case 14: // INTEGER
        value.move< int > (s.value);
        break;

      case 13: // IDENTIFIER
        value.move< std::string > (s.value);
        break;

      default:
        break;
    }

    location = s.location;
  }

  // by_type.
  inline
  OLAP_Parser::by_type::by_type ()
    : type (empty_symbol)
  {}

  inline
  OLAP_Parser::by_type::by_type (const by_type& other)
    : type (other.type)
  {}

  inline
  OLAP_Parser::by_type::by_type (token_type t)
    : type (yytranslate_ (t))
  {}

  inline
  void
  OLAP_Parser::by_type::clear ()
  {
    type = empty_symbol;
  }

  inline
  void
  OLAP_Parser::by_type::move (by_type& that)
  {
    type = that.type;
    that.clear ();
  }

  inline
  int
  OLAP_Parser::by_type::type_get () const
  {
    return type;
  }
  // Implementation of make_symbol for each symbol type.
  OLAP_Parser::symbol_type
  OLAP_Parser::make_BGN (const location_type& l)
  {
    return symbol_type (token::BGN, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_END (const location_type& l)
  {
    return symbol_type (token::END, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_CREATE (const location_type& l)
  {
    return symbol_type (token::CREATE, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_CUBE (const location_type& l)
  {
    return symbol_type (token::CUBE, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_LOAD (const location_type& l)
  {
    return symbol_type (token::LOAD, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_DROP (const location_type& l)
  {
    return symbol_type (token::DROP, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_COLUMN (const location_type& l)
  {
    return symbol_type (token::COLUMN, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_INFILE (const location_type& l)
  {
    return symbol_type (token::INFILE, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_AS (const location_type& l)
  {
    return symbol_type (token::AS, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_INTO (const location_type& l)
  {
    return symbol_type (token::INTO, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_IDENTIFIER (const std::string& v, const location_type& l)
  {
    return symbol_type (token::IDENTIFIER, v, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_INTEGER (const int& v, const location_type& l)
  {
    return symbol_type (token::INTEGER, v, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_HADAMARD (const location_type& l)
  {
    return symbol_type (token::HADAMARD, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_KRAO (const location_type& l)
  {
    return symbol_type (token::KRAO, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_KRON (const location_type& l)
  {
    return symbol_type (token::KRON, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_TR (const location_type& l)
  {
    return symbol_type (token::TR, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_VECTOR (const location_type& l)
  {
    return symbol_type (token::VECTOR, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_MATRIX (const location_type& l)
  {
    return symbol_type (token::MATRIX, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_BITMAP (const location_type& l)
  {
    return symbol_type (token::BITMAP, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_BANG (const location_type& l)
  {
    return symbol_type (token::BANG, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_TBL_READ (const location_type& l)
  {
    return symbol_type (token::TBL_READ, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_MX_FILTER_AND (const location_type& l)
  {
    return symbol_type (token::MX_FILTER_AND, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_TBL_WRITE (const location_type& l)
  {
    return symbol_type (token::TBL_WRITE, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_CONDITION (const location_type& l)
  {
    return symbol_type (token::CONDITION, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_KEY_CONDITION (const location_type& l)
  {
    return symbol_type (token::KEY_CONDITION, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_START (const location_type& l)
  {
    return symbol_type (token::START, l);
  }

  OLAP_Parser::symbol_type
  OLAP_Parser::make_STOP (const location_type& l)
  {
    return symbol_type (token::STOP, l);
  }



  // by_state.
  inline
  OLAP_Parser::by_state::by_state ()
    : state (empty_state)
  {}

  inline
  OLAP_Parser::by_state::by_state (const by_state& other)
    : state (other.state)
  {}

  inline
  void
  OLAP_Parser::by_state::clear ()
  {
    state = empty_state;
  }

  inline
  void
  OLAP_Parser::by_state::move (by_state& that)
  {
    state = that.state;
    that.clear ();
  }

  inline
  OLAP_Parser::by_state::by_state (state_type s)
    : state (s)
  {}

  inline
  OLAP_Parser::symbol_number_type
  OLAP_Parser::by_state::type_get () const
  {
    if (state == empty_state)
      return empty_symbol;
    else
      return yystos_[state];
  }

  inline
  OLAP_Parser::stack_symbol_type::stack_symbol_type ()
  {}


  inline
  OLAP_Parser::stack_symbol_type::stack_symbol_type (state_type s, symbol_type& that)
    : super_type (s, that.location)
  {
      switch (that.type_get ())
    {
      case 14: // INTEGER
        value.move< int > (that.value);
        break;

      case 13: // IDENTIFIER
        value.move< std::string > (that.value);
        break;

      default:
        break;
    }

    // that is emptied.
    that.type = empty_symbol;
  }

  inline
  OLAP_Parser::stack_symbol_type&
  OLAP_Parser::stack_symbol_type::operator= (const stack_symbol_type& that)
  {
    state = that.state;
      switch (that.type_get ())
    {
      case 14: // INTEGER
        value.copy< int > (that.value);
        break;

      case 13: // IDENTIFIER
        value.copy< std::string > (that.value);
        break;

      default:
        break;
    }

    location = that.location;
    return *this;
  }


  template <typename Base>
  inline
  void
  OLAP_Parser::yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yysym);
  }

#if YYDEBUG
  template <typename Base>
  void
  OLAP_Parser::yy_print_ (std::ostream& yyo,
                                     const basic_symbol<Base>& yysym) const
  {
    std::ostream& yyoutput = yyo;
    YYUSE (yyoutput);
    symbol_number_type yytype = yysym.type_get ();
    // Avoid a (spurious) G++ 4.8 warning about "array subscript is
    // below array bounds".
    if (yysym.empty ())
      std::abort ();
    yyo << (yytype < yyntokens_ ? "token" : "nterm")
        << ' ' << yytname_[yytype] << " ("
        << yysym.location << ": ";
    YYUSE (yytype);
    yyo << ')';
  }
#endif

  inline
  void
  OLAP_Parser::yypush_ (const char* m, state_type s, symbol_type& sym)
  {
    stack_symbol_type t (s, sym);
    yypush_ (m, t);
  }

  inline
  void
  OLAP_Parser::yypush_ (const char* m, stack_symbol_type& s)
  {
    if (m)
      YY_SYMBOL_PRINT (m, s);
    yystack_.push (s);
  }

  inline
  void
  OLAP_Parser::yypop_ (unsigned int n)
  {
    yystack_.pop (n);
  }

#if YYDEBUG
  std::ostream&
  OLAP_Parser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  OLAP_Parser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  OLAP_Parser::debug_level_type
  OLAP_Parser::debug_level () const
  {
    return yydebug_;
  }

  void
  OLAP_Parser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif // YYDEBUG

  inline OLAP_Parser::state_type
  OLAP_Parser::yy_lr_goto_state_ (state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - yyntokens_] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - yyntokens_];
  }

  inline bool
  OLAP_Parser::yy_pact_value_is_default_ (int yyvalue)
  {
    return yyvalue == yypact_ninf_;
  }

  inline bool
  OLAP_Parser::yy_table_value_is_error_ (int yyvalue)
  {
    return yyvalue == yytable_ninf_;
  }

  int
  OLAP_Parser::parse ()
  {
    // State.
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The locations where the error started and ended.
    stack_symbol_type yyerror_range[3];

    /// The return value of parse ().
    int yyresult;

    // FIXME: This shoud be completely indented.  It is not yet to
    // avoid gratuitous conflicts when merging into the master branch.
    try
      {
    YYCDEBUG << "Starting parse" << std::endl;


    /* Initialize the stack.  The initial state will be set in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystack_.clear ();
    yypush_ (YY_NULLPTR, 0, yyla);

    // A new symbol was pushed on the stack.
  yynewstate:
    YYCDEBUG << "Entering state " << yystack_[0].state << std::endl;

    // Accept?
    if (yystack_[0].state == yyfinal_)
      goto yyacceptlab;

    goto yybackup;

    // Backup.
  yybackup:

    // Try to take a decision without lookahead.
    yyn = yypact_[yystack_[0].state];
    if (yy_pact_value_is_default_ (yyn))
      goto yydefault;

    // Read a lookahead token.
    if (yyla.empty ())
      {
        YYCDEBUG << "Reading a token: ";
        try
          {
            yyla.type = yytranslate_ (yylex (&yyla.value, &yyla.location));
          }
        catch (const syntax_error& yyexc)
          {
            error (yyexc);
            goto yyerrlab1;
          }
      }
    YY_SYMBOL_PRINT ("Next token is", yyla);

    /* If the proper action on seeing token YYLA.TYPE is to reduce or
       to detect an error, take that action.  */
    yyn += yyla.type_get ();
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.type_get ())
      goto yydefault;

    // Reduce or error.
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
        if (yy_table_value_is_error_ (yyn))
          goto yyerrlab;
        yyn = -yyn;
        goto yyreduce;
      }

    // Count tokens shifted since error; after three, turn off error status.
    if (yyerrstatus_)
      --yyerrstatus_;

    // Shift the lookahead token.
    yypush_ ("Shifting", yyn, yyla);
    goto yynewstate;

  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[yystack_[0].state];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;

  /*-----------------------------.
  | yyreduce -- Do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    {
      stack_symbol_type yylhs;
      yylhs.state = yy_lr_goto_state_(yystack_[yylen].state, yyr1_[yyn]);
      /* Variants are always initialized to an empty instance of the
         correct type. The default '$$ = $1' action is NOT applied
         when using variants.  */
        switch (yyr1_[yyn])
    {
      case 14: // INTEGER
        yylhs.value.build< int > ();
        break;

      case 13: // IDENTIFIER
        yylhs.value.build< std::string > ();
        break;

      default:
        break;
    }


      // Compute the default @$.
      {
        slice<stack_symbol_type, stack_type> slice (yystack_, yylen);
        YYLLOC_DEFAULT (yylhs.location, slice, yylen);
      }

      // Perform the reduction.
      YY_REDUCE_PRINT (yyn);
      try
        {
          switch (yyn)
            {
  case 10:
#line 80 "olap_parser.yy" // lalr1.cc:859
    {
                   std::cout << "created cube " << yystack_[0].value.as< std::string > () << std::endl;
}
#line 870 "olap_parser.tab.cc" // lalr1.cc:859
    break;

  case 11:
#line 85 "olap_parser.yy" // lalr1.cc:859
    {
driver.load_matrix_csc( yystack_[4].value.as< std::string > (), yystack_[6].value.as< int > ());
}
#line 878 "olap_parser.tab.cc" // lalr1.cc:859
    break;

  case 12:
#line 88 "olap_parser.yy" // lalr1.cc:859
    {
 }
#line 885 "olap_parser.tab.cc" // lalr1.cc:859
    break;

  case 13:
#line 93 "olap_parser.yy" // lalr1.cc:859
    {
}
#line 892 "olap_parser.tab.cc" // lalr1.cc:859
    break;

  case 14:
#line 96 "olap_parser.yy" // lalr1.cc:859
    {
}
#line 899 "olap_parser.tab.cc" // lalr1.cc:859
    break;

  case 18:
#line 105 "olap_parser.yy" // lalr1.cc:859
    {
     }
#line 906 "olap_parser.tab.cc" // lalr1.cc:859
    break;

  case 19:
#line 107 "olap_parser.yy" // lalr1.cc:859
    {

}
#line 914 "olap_parser.tab.cc" // lalr1.cc:859
    break;


#line 918 "olap_parser.tab.cc" // lalr1.cc:859
            default:
              break;
            }
        }
      catch (const syntax_error& yyexc)
        {
          error (yyexc);
          YYERROR;
        }
      YY_SYMBOL_PRINT ("-> $$ =", yylhs);
      yypop_ (yylen);
      yylen = 0;
      YY_STACK_PRINT ();

      // Shift the result of the reduction.
      yypush_ (YY_NULLPTR, yylhs);
    }
    goto yynewstate;

  /*--------------------------------------.
  | yyerrlab -- here on detecting error.  |
  `--------------------------------------*/
  yyerrlab:
    // If not already recovering from an error, report this error.
    if (!yyerrstatus_)
      {
        ++yynerrs_;
        error (yyla.location, yysyntax_error_ (yystack_[0].state, yyla));
      }


    yyerror_range[1].location = yyla.location;
    if (yyerrstatus_ == 3)
      {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.type_get () == yyeof_)
          YYABORT;
        else if (!yyla.empty ())
          {
            yy_destroy_ ("Error: discarding", yyla);
            yyla.clear ();
          }
      }

    // Else will try to reuse lookahead token after shifting the error token.
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:

    /* Pacify compilers like GCC when the user code never invokes
       YYERROR and the label yyerrorlab therefore never appears in user
       code.  */
    if (false)
      goto yyerrorlab;
    yyerror_range[1].location = yystack_[yylen - 1].location;
    /* Do not reclaim the symbols of the rule whose action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    goto yyerrlab1;

  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;   // Each real token shifted decrements this.
    {
      stack_symbol_type error_token;
      for (;;)
        {
          yyn = yypact_[yystack_[0].state];
          if (!yy_pact_value_is_default_ (yyn))
            {
              yyn += yyterror_;
              if (0 <= yyn && yyn <= yylast_ && yycheck_[yyn] == yyterror_)
                {
                  yyn = yytable_[yyn];
                  if (0 < yyn)
                    break;
                }
            }

          // Pop the current state because it cannot handle the error token.
          if (yystack_.size () == 1)
            YYABORT;

          yyerror_range[1].location = yystack_[0].location;
          yy_destroy_ ("Error: popping", yystack_[0]);
          yypop_ ();
          YY_STACK_PRINT ();
        }

      yyerror_range[2].location = yyla.location;
      YYLLOC_DEFAULT (error_token.location, yyerror_range, 2);

      // Shift the error token.
      error_token.state = yyn;
      yypush_ ("Shifting", error_token);
    }
    goto yynewstate;

    // Accept.
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;

    // Abort.
  yyabortlab:
    yyresult = 1;
    goto yyreturn;

  yyreturn:
    if (!yyla.empty ())
      yy_destroy_ ("Cleanup: discarding lookahead", yyla);

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    while (1 < yystack_.size ())
      {
        yy_destroy_ ("Cleanup: popping", yystack_[0]);
        yypop_ ();
      }

    return yyresult;
  }
    catch (...)
      {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack"
                 << std::endl;
        // Do not try to display the values of the reclaimed symbols,
        // as their printer might throw an exception.
        if (!yyla.empty ())
          yy_destroy_ (YY_NULLPTR, yyla);

        while (1 < yystack_.size ())
          {
            yy_destroy_ (YY_NULLPTR, yystack_[0]);
            yypop_ ();
          }
        throw;
      }
  }

  void
  OLAP_Parser::error (const syntax_error& yyexc)
  {
    error (yyexc.location, yyexc.what());
  }

  // Generate an error message.
  std::string
  OLAP_Parser::yysyntax_error_ (state_type, const symbol_type&) const
  {
    return YY_("syntax error");
  }


  const signed char OLAP_Parser::yypact_ninf_ = -14;

  const signed char OLAP_Parser::yytable_ninf_ = -1;

  const signed char
  OLAP_Parser::yypact_[] =
  {
       2,    -1,    13,    15,   -13,    10,    18,    19,   -14,   -14,
      -4,   -14,     3,     4,     5,    23,   -14,    24,    29,    30,
     -14,   -14,   -14,   -14,   -14,   -14,   -14,   -14,     9,     1,
      11,   -14,    28,    31,   -11,   -14,    14,   -14,    33,    36,
      -6,   -11,   -14,   -14,    34,    35,    37,    38,    39,   -14,
      40,    20,    44,    45,   -14,   -14,   -14,   -14,   -14,    46,
      47,    49,    50,    51,    52,   -14,   -14
  };

  const unsigned char
  OLAP_Parser::yydefact_[] =
  {
       0,     5,     0,     0,     0,     0,     0,     0,    18,    19,
       0,     4,     0,     0,     0,     0,     1,     0,     0,     0,
      13,    14,    15,     2,     3,     6,     7,     8,     0,     0,
       0,    10,     0,     0,     0,     9,     0,    17,     0,     0,
       0,     0,    20,    16,     0,     0,     0,     0,     0,    25,
       0,     0,     0,     0,    22,    24,    23,    21,    26,     0,
       0,     0,     0,     0,     0,    11,    12
  };

  const signed char
  OLAP_Parser::yypgoto_[] =
  {
     -14,   -14,   -14,    48,   -14,   -14,   -14,   -14,    41,    42,
       8
  };

  const signed char
  OLAP_Parser::yydefgoto_[] =
  {
      -1,     2,    10,    11,    12,    13,    14,    29,    15,    30,
      42
  };

  const unsigned char
  OLAP_Parser::yytable_[] =
  {
      23,     3,    40,     4,     3,     1,     4,    18,    19,    46,
      47,    48,    49,    16,    28,     5,     6,     7,     5,     6,
       7,    17,    41,    20,     8,     9,    50,     8,     9,     8,
       9,    21,    22,    25,    26,    27,    28,    31,    32,    33,
      34,    37,    38,    44,    43,    39,    45,    52,    53,    51,
      54,    55,    56,    57,    58,    59,    60,     0,    24,    61,
      62,    63,    64,     0,    65,    66,     0,     0,     0,     0,
      35,    36
  };

  const signed char
  OLAP_Parser::yycheck_[] =
  {
       4,     5,    13,     7,     5,     3,     7,    20,    21,    15,
      16,    17,    18,     0,    13,    19,    20,    21,    19,    20,
      21,     6,    33,    13,    28,    29,    32,    28,    29,    28,
      29,    13,    13,    30,    30,    30,    13,    13,     9,     9,
      31,    30,    14,    10,    30,    14,    10,    13,    13,    41,
      13,    13,    13,    13,    34,    11,    11,    -1,    10,    13,
      13,    12,    12,    -1,    13,    13,    -1,    -1,    -1,    -1,
      29,    29
  };

  const unsigned char
  OLAP_Parser::yystos_[] =
  {
       0,     3,    36,     5,     7,    19,    20,    21,    28,    29,
      37,    38,    39,    40,    41,    43,     0,     6,    20,    21,
      13,    13,    13,     4,    38,    30,    30,    30,    13,    42,
      44,    13,     9,     9,    31,    43,    44,    30,    14,    14,
      13,    33,    45,    30,    10,    10,    15,    16,    17,    18,
      32,    45,    13,    13,    13,    13,    13,    13,    34,    11,
      11,    13,    13,    12,    12,    13,    13
  };

  const unsigned char
  OLAP_Parser::yyr1_[] =
  {
       0,    35,    36,    37,    37,    37,    38,    38,    38,    38,
      39,    40,    40,    41,    41,    41,    42,    42,    43,    43,
      44,    45,    45,    45,    45,    45,    45
  };

  const unsigned char
  OLAP_Parser::yyr2_[] =
  {
       0,     2,     3,     2,     1,     0,     2,     2,     2,     3,
       3,    10,    10,     2,     2,     2,     3,     2,     1,     1,
       3,     3,     3,     3,     3,     2,     3
  };


#if YYDEBUG
  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a yyntokens_, nonterminals.
  const char*
  const OLAP_Parser::yytname_[] =
  {
  "$end", "error", "$undefined", "BGN", "END", "CREATE", "CUBE", "LOAD",
  "DROP", "COLUMN", "INFILE", "AS", "INTO", "IDENTIFIER", "INTEGER",
  "HADAMARD", "KRAO", "KRON", "TR", "VECTOR", "MATRIX", "BITMAP", "BANG",
  "TBL_READ", "MX_FILTER_AND", "TBL_WRITE", "CONDITION", "KEY_CONDITION",
  "START", "STOP", "';'", "'='", "'*'", "'('", "')'", "$accept",
  "initial_expression", "body", "elem", "Create_declaration",
  "Load_declaration", "Inline_declaration_type", "query", "time",
  "atribuition_expression", "expression", YY_NULLPTR
  };


  const unsigned char
  OLAP_Parser::yyrline_[] =
  {
       0,    66,    66,    69,    70,    71,    74,    75,    76,    77,
      80,    85,    88,    92,    95,    98,   101,   102,   105,   107,
     113,   116,   117,   118,   119,   120,   121
  };

  // Print the state stack on the debug stream.
  void
  OLAP_Parser::yystack_print_ ()
  {
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator
           i = yystack_.begin (),
           i_end = yystack_.end ();
         i != i_end; ++i)
      *yycdebug_ << ' ' << i->state;
    *yycdebug_ << std::endl;
  }

  // Report on the debug stream that the rule \a yyrule is going to be reduced.
  void
  OLAP_Parser::yy_reduce_print_ (int yyrule)
  {
    unsigned int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
               << " (line " << yylno << "):" << std::endl;
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
                       yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // YYDEBUG

  // Symbol number corresponding to token number t.
  inline
  OLAP_Parser::token_number_type
  OLAP_Parser::yytranslate_ (int t)
  {
    static
    const token_number_type
    translate_table[] =
    {
     0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      33,    34,    32,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,    30,
       2,    31,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29
    };
    const unsigned int user_token_number_max_ = 284;
    const token_number_type undef_token_ = 2;

    if (static_cast<int>(t) <= yyeof_)
      return yyeof_;
    else if (static_cast<unsigned int> (t) <= user_token_number_max_)
      return translate_table[t];
    else
      return undef_token_;
  }

#line 12 "olap_parser.yy" // lalr1.cc:1167
} // OLAP
#line 1287 "olap_parser.tab.cc" // lalr1.cc:1167
#line 124 "olap_parser.yy" // lalr1.cc:1168



void
OLAP::OLAP_Parser::error( const location_type &l, const std::string &err_message )
    {
std::cerr << "Error: " << err_message << " at " << l << "\n";
   }

