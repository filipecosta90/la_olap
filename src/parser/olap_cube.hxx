#ifndef __OLAP_CUBE_HXX__
#define __OLAP_CUBE_HXX__ 1

namespace OLAP{
  class OLAP_Cube {
    private:

    public:
      OLAP_Cube( std::string cube_name ); // constructor
      int get_row_from_string(std::string field );
      bool load_matrix_to_csc_from_tbl( std::string filename, int col_number, int max_col_size );
  };
}/* end namespace OLAP */
#endif /* END __OLAP_CUBE_HXX__ */
