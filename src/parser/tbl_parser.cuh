#ifndef __TBL_PARSER_CUH__
#define __TBL_PARSER_CUH__ 1


namespace TBL{
  class TBL_Parser {
    void load_matrix_csc ( 
        std::string filename, int col_number, int max_col_size, 
        int* n_nnz, int* n_rows, int* n_cols,
        float** __restrict__  A_csc_values,
        int** __restrict__  A_row_ind,
        int** __restrict__  A_col_ptr
        );
  };

}/* end namespace TBL */
#endif /* END __TBL_PARSER_CUH__ */
