#include <stdio.h>
#include <glib.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mkl_spblas.h"
#include "mkl_types.h"
#include "mkl.h"
#include "olap.h"

//Cache-Lines size is (typically) 64 bytes
#define MEM_LINE_SIZE 64
#define ARRAY_SIZE MEM_LINE_SIZE / sizeof (MKL_INT)
#define GROWTH_FACTOR 2
#define MAX_FIELD_SIZE 128
#define MAX_REG_SIZE 1024

int main( int argc, char* argv[]){
  mkl_verbose(1);

  //define COO sparse-matrix M
  float* A_csr_values = NULL;
  float* B_csr_values = NULL;
  float* C_csr_values = NULL;

  MKL_INT* A_JA;
  MKL_INT* A_IA;
  MKL_INT* B_JA;
  MKL_INT* B_IA;
  MKL_INT* C_JA;
  MKL_INT* C_IA;

  MKL_INT A_rows;
  MKL_INT A_columns;
  MKL_INT A_nnz;

  MKL_INT B_rows;
  MKL_INT B_columns;
  MKL_INT B_nnz;

  MKL_INT C_rows;
  MKL_INT C_columns;
  MKL_INT C_nnz;

  MKL_INT tbl_column = atoi (argv[2]);
  MKL_INT tbl_column_2 = atoi (argv[3]);
  //read A
  tbl_read( argv[1], tbl_column, &A_nnz, &A_rows, &A_columns , &A_csr_values, &A_JA, &A_IA);
  print_csr( A_csr_values, A_JA, A_IA, A_nnz, A_rows, A_columns);

  //read B
  tbl_read( argv[1], tbl_column_2, &B_nnz, &B_rows, &B_columns , &B_csr_values, &B_JA, &B_IA);
  print_csr( B_csr_values, B_JA, B_IA, B_nnz, B_rows, B_columns);

  //   compute C = A hadarmard B
  //  csr_hadamard( A_csr_values, A_JA, A_IA, A_nnz, A_rows, B_csr_values, B_JA, B_IA , B_nnz, &C_csr_values, &C_JA, &C_IA, &C_nnz );
  //  print_csr( C_csr_values, C_JA, C_IA, C_nnz, B_rows, B_columns);

  // compute C = A krao B
  csr_krao( 
      A_csr_values, A_JA, A_IA, A_nnz, A_rows, A_columns, 
      B_csr_values, B_JA, B_IA, B_nnz, A_rows, B_columns, 
      &C_csr_values, &C_JA, &C_IA, &C_nnz, &C_rows, &C_columns 
      );
  print_csr( C_csr_values, C_JA, C_IA, C_nnz, C_rows, C_columns);

  /*
  /////////////////////////////////
  //
  //   CONVERT FROM CSR TO BSR
  //
  ////////////////////////////////

  sparse_matrix_t  A_bsr;
  sparse_status_t status_convert_bsr;
  MKL_INT mblk = 1;

  printf("going to convert to BSR:\n");
  status_convert_bsr = mkl_sparse_convert_bsr ( A_csr,   mblk, SPARSE_LAYOUT_ROW_MAJOR, SPARSE_OPERATION_NON_TRANSPOSE, &A_bsr );
  check_errors(status_convert_bsr); 

  float* bsr_values = NULL;
  MKL_INT* JA_B;
  MKL_INT* IA_B;
  MKL_INT* ROWS_END;
  MKL_INT ROWS_B = 0;
  MKL_INT COLS_B = 0;
  MKL_INT BLOCK_SIZE;

  sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO; 
  sparse_layout_t block_layout = SPARSE_LAYOUT_ROW_MAJOR;

  sparse_status_t status_export_bsr;
  printf("going to export BSR:\n");
  mkl_verbose(1);
  status_export_bsr = mkl_sparse_s_export_bsr ( A_bsr, &indexing, &block_layout, &ROWS_B, &COLS_B, &BLOCK_SIZE, &IA_B, &ROWS_END, &JA_B, &bsr_values);
  check_errors(status_export_bsr); 
  mkl_verbose(1);
  for ( int i=0; i < NNZ; i++){
  printf("%f, ", bsr_values[i]);
  }

  printf("\n");
  for ( int i=0; i < mblk; i++){ 
  printf("%d, ", JA_B[i]);
  }
  printf("\n");
  for ( int i=0; i < ROWS_B+1; i++){
  printf("%d, ", IA_B[i]);
  }
  printf("\n");
  */
  return 0;
}
