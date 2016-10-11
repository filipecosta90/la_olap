#include <cctype>
#include <fstream>
#include <cassert>

#include "olap_driver.hpp"

OLAP::OLAP_Driver::~OLAP_Driver()
{
    delete(scanner);
    scanner = nullptr;
    delete(parser);
    parser = nullptr;
}

void
OLAP::OLAP_Driver::parse( const char * const filename )
{
    assert( filename != nullptr );
    std::ifstream in_file( filename );
    if( ! in_file.good() )
    {
        exit( EXIT_FAILURE );
    }
    parse_helper( in_file );
    return;
}

void
OLAP::OLAP_Driver::parse( std::istream &stream )
{
    if( ! stream.good()  && stream.eof() )
    {
        return;
    }
    //else
    parse_helper( stream );
    return;
}


void
OLAP::OLAP_Driver::parse_helper( std::istream &stream )
{
    
    delete(scanner);
    try
    {
        scanner = new OLAP::OLAP_Scanner( &stream );
    }
    catch( std::bad_alloc &ba )
    {
        std::cerr << "Failed to allocate scanner: (" <<
        ba.what() << "), exiting!!\n";
        exit( EXIT_FAILURE );
    }
    
    delete(parser);
    try
    {
        parser = new OLAP::OLAP_Parser( (*scanner) /* scanner */,
                                   (*this) /* driver */ );
    }
    catch( std::bad_alloc &ba )
    {
        std::cerr << "Failed to allocate parser: (" <<
        ba.what() << "), exiting!!\n";
        exit( EXIT_FAILURE );
    }
    const int accept( 0 );
    if( parser->parse() != accept )
    {
        std::cerr << "Parse failed!!\n";
    }
    return;
}

void
OLAP::OLAP_Driver::col_read_csc ( std::string  filename, int col_number ){
}


std::ostream&
OLAP::OLAP_Driver::print( std::ostream &stream )
{
    stream << red  << "Debug info: " << norm << "\n";
    stream << blue << "OLAP: "  << norm << "\n";
    
    return(stream);
}
