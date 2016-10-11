#ifndef __OLAPDRIVER_HPP__
#define __OLAPDRIVER_HPP__ 1

#include <string>
#include <cstddef>
#include <istream>

#include "olap_scanner.hpp"
#include "olap_parser.tab.hh"

namespace OLAP{
    
    class OLAP_Driver{
    public:
        OLAP_Driver() = default;
        
        virtual ~OLAP_Driver();
        
        /**
         * parse - parse from a file
         * @param filename - valid string with input file
         */
        void parse( const char * const filename );
        /**
         * parse - parse from a c++ input stream
         * @param is - std::istream&, valid input stream
         */
        void parse( std::istream &iss );
        
        void add_upper();
        void add_lower();
        void add_word( const std::string &word );
        void add_newline();
        void add_char();
        
        std::ostream& print(std::ostream &stream);
    private:
        
        void parse_helper( std::istream &stream );
        
        std::size_t  chars      = 0;
        std::size_t  words      = 0;
        std::size_t  lines      = 0;
        std::size_t  uppercase  = 0;
        std::size_t  lowercase  = 0;
        OLAP::OLAP_Parser  *parser  = nullptr;
        OLAP::OLAP_Scanner *scanner = nullptr;
        
        const std::string red   = "\033[1;31m";
        const std::string blue  = "\033[1;36m";
        const std::string norm  = "\033[0m";
    };
    
} /* end namespace OLAP */
#endif /* END __OLAPDRIVER_HPP__ */
