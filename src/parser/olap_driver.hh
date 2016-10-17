#ifndef __OLAPDRIVER_HH__
#define __OLAPDRIVER_HH__ 1

#include <string>
#include <cstddef>
#include <istream>

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

      void load_matrix_csc ( std::string  filename, int col_number );

      std::ostream& print(std::ostream &stream);
    private:

      void parse_helper( std::istream &stream );

      OLAP::OLAP_Parser  *parser  = nullptr;
      OLAP::OLAP_Scanner *scanner = nullptr;

      const std::string red   = "\033[1;31m";
      const std::string blue  = "\033[1;36m";
      const std::string norm  = "\033[0m";
  };

} /* end namespace OLAP */
#endif /* END __OLAPDRIVER_HPP__ */
