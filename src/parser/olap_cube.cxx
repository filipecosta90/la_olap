#include <cctype>
#include <fstream>
#include <cassert>
#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <string>

//GLIB
#include <glib.h>

#include "olap_cube.hxx"

OLAP::OLAP_Cube::OLAP_Cube( std::string cube_name ){

}

int OLAP::OLAP_Cube::get_row_from_string(std::string field ){
  // since quarks start by 1 and we want to start at line 0 and not 1 lets decrement
  return (((int) g_quark_from_string ( field.c_str() )) -1 );
}

bool OLAP::OLAP_Cube::load_matrix_to_csc_from_tbl( std::string filename, int col_number, int max_col_size ){
  /* 
     int n_nnz, n_rows, n_cols;
     float *A_csc_values;
     int *A_row_ind, *A_col_ptr;
     driver.load_matrix_csc( $6, $4, 100, &n_nnz, &n_rows, &n_cols, &A_csc_values, &A_row_ind, &A_col_ptr );
     std::cout << "NNZ " << n_nnz << "|\tN_ROWS " << n_rows << "|\tN_COLS " << n_cols << std::endl;
     */
  return true;
}

