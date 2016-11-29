#ifndef __OLAPDRIVER_HH__
#define __OLAPDRIVER_HH__ 1

#include <string>
#include <iostream>
#include <fstream>
#include <cstddef>

#define MEM_LINE_SIZE 32

namespace OLAP{
  // Forward declarations of classes
  class OLAP_Parser;
  class OLAP_Scanner;
  class OLAP_Engine;
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

      void load_matrix_csc ( 
          std::string  filename, int col_number, int max_col_size,
          int* n_nnz, int* n_rows, int* n_cols,
          float** __restrict__ A_csc_values,
          int** __restrict__ A_row_ind,
          int** __restrict__ A_col_ptr
          );

      std::ostream& print(std::ostream &stream);
    private:

      void parse_helper( std::istream &stream );

      OLAP::OLAP_Parser  *parser  = nullptr;
      OLAP::OLAP_Scanner *scanner = nullptr;
      OLAP::OLAP_Engine *engine = nullptr;
      /// Allows Parser and Scanner to access private attributes
      //            /// of the Driver class
      friend class  OLAP_Parser;
      friend class  OLAP_Scanner;
      friend class  OLAP_Engine;

      const std::string red   = "\033[1;31m";
      const std::string blue  = "\033[1;36m";
      const std::string norm  = "\033[0m";
  };

} /* end namespace OLAP */
#endif /* END __OLAPDRIVER_HH__ */
