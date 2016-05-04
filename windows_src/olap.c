#include <stdio.h>
#include <glib.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mkl_spblas.h"
#include "mkl_types.h"
#include "mkl.h"


//Cache-Lines size is (typically) 64 bytes
#define MEM_LINE_SIZE 64
#define ARRAY_SIZE MEM_LINE_SIZE / sizeof (MKL_INT)
#define GROWTH_FACTOR 2
#define MAX_FIELD_SIZE 128
#define MAX_REG_SIZE 1024

char* getfield(char* line, int num, char* return_string ){
  return_string = strtok(line, "|\n");
  int pos = 1;
  for ( ; pos < num; pos++ ){
    return_string = strtok(NULL, "|\n");
  }
  return return_string;
}

void check_errors( sparse_status_t stat ){
  if ( stat == SPARSE_STATUS_SUCCESS ){
    printf( "SPARSE_STATUS_SUCCESS.\n");
  }
  if ( stat == SPARSE_STATUS_NOT_INITIALIZED ){
    printf( "SPARSE_STATUS_NOT_INITIALIZED.\n");
  }
  if ( stat == SPARSE_STATUS_ALLOC_FAILED ){
    printf( "SPARSE_STATUS_ALLOC_FAILED.\n");
  }
  if ( stat == SPARSE_STATUS_INVALID_VALUE ){
    printf( "SPARSE_STATUS_INVALID_VALUE.\n");
  }
  if ( stat == SPARSE_STATUS_EXECUTION_FAILED){
    printf( "SPARSE_STATUS_EXECUTION_FAILED.\n");
  }
  if ( stat == SPARSE_STATUS_INTERNAL_ERROR){
    printf( "SPARSE_STATUS_INTERNAL_ERROR.\n");
  }
  if ( stat == SPARSE_STATUS_NOT_SUPPORTED){
    printf( "SPARSE_STATUS_NOT_SUPPORTED.\n");
  }
}

int main( int argc, char* argv[]){
  mkl_verbose(1);
  MKL_INT current_values_size = ARRAY_SIZE;

  //define COO sparse-matrix M
  MKL_INT* aux_coo_rows;
  aux_coo_rows = (MKL_INT*) malloc (current_values_size * sizeof(MKL_INT));

  FILE* stream = fopen(argv[1], "r");
  char line[1024];
  MKL_INT column = atoi (argv[2]);
  MKL_INT number_rows = - 1;
  MKL_INT number_columns = -1 ;
  MKL_INT element_number = 1;
  MKL_INT job[8];

  for( element_number = 0 ; (fgets(line, MAX_REG_SIZE, stream) ) ; ++element_number )
  {
    char *field = (char*) malloc( MAX_FIELD_SIZE * sizeof(char) );
    char* tmp_field = strdup(line);
    field = getfield(tmp_field, column, field);
    MKL_INT quark_field;

    quark_field = g_quark_from_string (field);
    printf("%s %d\n", field, quark_field-1);
    /* if arrays are full double its size */
    if ( element_number >= current_values_size ){
      current_values_size *= GROWTH_FACTOR;
      aux_coo_rows = (MKL_INT*) realloc(aux_coo_rows, (current_values_size) * GROWTH_FACTOR * sizeof(MKL_INT) );
    }

    /* normal coo property */
    aux_coo_rows[element_number]=  quark_field -1 ;

    if (  quark_field > number_rows ) {
      number_rows =  quark_field;
    }
    free(tmp_field);
  }

  MKL_INT NNZ = element_number;
  number_columns = element_number;

  //define COO sparse-matrix M
  float* coo_values;
  MKL_INT* coo_rows;
  MKL_INT* coo_columns;

  coo_values = (float*) mkl_malloc ((NNZ * sizeof(float)), MEM_LINE_SIZE );
  coo_rows =  (MKL_INT*) mkl_malloc ((NNZ * sizeof(MKL_INT)), MEM_LINE_SIZE );
  coo_columns =  (MKL_INT*) mkl_malloc ((NNZ * sizeof(MKL_INT)), MEM_LINE_SIZE );

  for (int pos = 0; pos < NNZ; pos++) {
    coo_values[pos] = 1.0;
    coo_columns[pos] = pos;
    coo_rows[pos] = aux_coo_rows[pos];
    printf("%f %d %d\n", coo_values[pos], coo_rows[pos], coo_columns[pos]);
  }
  free(aux_coo_rows);

  /////////////////////////////////
  //
  //   CONVERT FROM COO TO CSR
  //
  ////////////////////////////////

  float* csr_values = NULL;

  job[0]=  2; // if job[0]=2, the matrix in the coordinate format is converted to the CSR
  // format, and the column indices in CSR representation are sorted in the
  // increasing order within each row.
  job[1]= 0; // If job[1]=0, zero-based indexing for the matrix in CSR format is used;
  job[2]= 0; // If job[2]=0, zero-based indexing for the matrix in coordinate format is
  //used;
  job[3]= 0;
  job[4]= NNZ;// job[4]=nzmax - maximum number of the non-zero elements allowed if
  // job[0]=0.
  job[5]= 0;   // If job[5]=0, all arrays acsr, ja, ia are filled in for the output storage.

  MKL_INT* JA;
  MKL_INT* IA;

  csr_values = (float*) mkl_malloc ((element_number * sizeof(float)), MEM_LINE_SIZE );
  JA = (MKL_INT*) mkl_malloc (( element_number * sizeof(MKL_INT)), MEM_LINE_SIZE );
  IA = (MKL_INT*) mkl_malloc ((number_rows+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );
  printf("going to convert matrix N x M = %d x %d\nNNZ: %d JOB4: %d\n", number_rows, number_columns, NNZ, job[4]);
  sparse_status_t status_coo_csr;
  mkl_scsrcoo (job, &number_rows, csr_values, JA, IA, &NNZ, coo_values, coo_rows, coo_columns, &status_coo_csr);
  check_errors(status_coo_csr);

  for (MKL_INT pos = 0; pos < NNZ; pos++){
    printf("%f, ", csr_values[pos]);
  }
  printf("\n");
  for (int pos = 0; pos < NNZ; pos++){
    printf("%d, ", JA[pos]);
  }
  printf("\n");
  for (int pos = 0; pos <= number_rows; pos++){
    printf("%d, ", IA[pos]);
  }
  printf("\n");

  mkl_free(coo_values);
  mkl_free(coo_rows);
  mkl_free(coo_columns);


  /////////////////////////////////
  //
  //   CREATE CSR
  //
  ////////////////////////////////

  sparse_matrix_t  A_csr;
  sparse_status_t status_create_csr;

  printf("going to create CSR:\n");
  status_create_csr = mkl_sparse_s_create_csr ( &A_csr, SPARSE_INDEX_BASE_ZERO, number_rows, number_columns, IA, IA+1, JA, csr_values);
  check_errors(status_create_csr); 

  /////////////////////////////////
  //
  //   TRANSPOSE CSR
  //
  ////////////////////////////////

  // If job[0]=0, the matrix in the CSR format is converted to the CSC format;
  job[0] = 0;
  // job[1]
  job[1] = 0;
  // If job[1]=0, zero-based indexing for the matrix in CSR format is used;
  // if job[1]=1, one-based indexing for the matrix in CSR format is used.
  // job[2]
  // If job[2]=0, zero-based indexing for the matrix in the CSC format is used;
  // if job[2]=1, one-based indexing for the matrix in the CSC format is used.
  job[2] = 0;
  // job[5] - job indicator.
  // If job[5]=0, only arrays ja1, ia1 are filled in for the output storage.
  // If job[5]≠0, all output arrays acsc, ja1, and ia1 are filled in for the output storage.
  job[5] = 1;


  float* csc_values = NULL;
  MKL_INT* JA1;
  MKL_INT* IA1;

  csc_values = (float*) mkl_malloc ((element_number * sizeof(float)), MEM_LINE_SIZE );
  JA1 = (MKL_INT*) mkl_malloc (( element_number * sizeof(MKL_INT)), MEM_LINE_SIZE );
  IA1 = (MKL_INT*) mkl_malloc ((number_columns+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );
  sparse_status_t status_convert_csc;

  printf("going to transpose CSR:\n");
  mkl_scsrcsc(job, &NNZ, csr_values, JA, IA, csc_values, JA1, IA1, &status_convert_csc);
  check_errors(status_convert_csc); 

  for (MKL_INT pos = 0; pos < NNZ; pos++){
    printf("%f, ", csc_values[pos]);
  }
  printf("\n");
  for (int pos = 0; pos < NNZ; pos++){
    printf("%d, ", JA1[pos]);
  }
  printf("\n");
  for (int pos = 0; pos <= number_columns; pos++){
    printf("%d, ", IA1[pos]);
  }
  printf("\n");

  /////////////////////////////////
  //
  //   Khatri-Rao CSR 
  //   C = A krao B 
  //
  ////////////////////////////////

  // If job[0]=0, the matrix in the CSR format is converted to the CSC format;
  job[0] = 0;
  // job[1]
  job[1] = 0;
  // If job[1]=0, zero-based indexing for the matrix in CSR format is used;
  // if job[1]=1, one-based indexing for the matrix in CSR format is used.
  // job[2]
  // If job[2]=0, zero-based indexing for the matrix in the CSC format is used;
  // if job[2]=1, one-based indexing for the matrix in the CSC format is used.
  job[2] = 0;
  // job[5] - job indicator.
  // If job[5]=0, only arrays ja1, ia1 are filled in for the output storage.
  // If job[5]≠0, all output arrays acsc, ja1, and ia1 are filled in for the output storage.
  job[5] = 1;

  float* csc_values_A = NULL;
  MKL_INT* JA1_A;
  MKL_INT* IA1_A;

  csc_values_A = (float*) mkl_malloc ((element_number * sizeof(float)), MEM_LINE_SIZE );
  JA1_A = (MKL_INT*) mkl_malloc (( element_number * sizeof(MKL_INT)), MEM_LINE_SIZE );
  IA1_A = (MKL_INT*) mkl_malloc ((number_columns+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );

  float* csc_values_B = NULL;
  MKL_INT* JA1_B;
  MKL_INT* IA1_B;

  csc_values_B = (float*) mkl_malloc ((element_number * sizeof(float)), MEM_LINE_SIZE );
  JA1_B = (MKL_INT*) mkl_malloc (( element_number * sizeof(MKL_INT)), MEM_LINE_SIZE );
  IA1_B = (MKL_INT*) mkl_malloc ((number_columns+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );

  printf("going to convert CSR A to CSC A:\n");
  mkl_scsrcsc(job, &number_columns, csr_values, JA, IA, csc_values_A, JA1_A, IA1_A, &status_convert_csc);
  check_errors(status_convert_csc); 

  printf("going to convert CSR B to CSC B:\n");
  mkl_scsrcsc(job, &number_columns, csr_values, JA, IA, csc_values_B, JA1_B, IA1_B, &status_convert_csc);
  check_errors(status_convert_csc); 

  MKL_INT number_rows_A = number_rows;
  MKL_INT number_cols_A = number_columns;
  MKL_INT number_cols_C = number_columns;
  MKL_INT number_rows_B = number_rows;

  MKL_INT number_rows_C = number_rows_A * number_rows_B;
  MKL_INT nnz_A = NNZ;
  MKL_INT nnz_B = NNZ;
  MKL_INT nnz_C = nnz_B * number_rows_A;

  float* csc_values_C = NULL;
  MKL_INT* JA1_C;
  MKL_INT* IA1_C;

  csc_values_C = (float*) mkl_malloc (( nnz_C * sizeof(float)), MEM_LINE_SIZE );
  JA1_C = (MKL_INT*) mkl_malloc (( nnz_C * sizeof(MKL_INT)), MEM_LINE_SIZE );
  IA1_C = (MKL_INT*) mkl_malloc ((number_cols_C+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );

  // go through all columns of B 
  // lets refer to it as column X of matrix B
  MKL_INT lower_limit_A, upper_limit_A, lower_limit_B, upper_limit_B; 
  MKL_INT lower_pos_A, upper_pos_A, lower_pos_B, upper_pos_B; 
  MKL_INT c_value_pos = 0;
  for (int a_column_pos = 0; a_column_pos < number_cols_A ; a_column_pos++ ){
    lower_pos_A = JA1_A[a_column_pos];
    upper_pos_A = JA1_A[a_column_pos+1] ;
    lower_pos_B = JA1_B[a_column_pos];
    upper_pos_B = JA1_B[a_column_pos+1] ;
    lower_limit_A = IA1_A[lower_pos_A];
    upper_limit_A = IA1_A[upper_pos_A];
    lower_limit_B = IA1_B[lower_pos_B];
    upper_limit_B = IA1_B[upper_pos_B];
    printf("col %d :: A( %d:%d) B(%d:%d)\n", a_column_pos, lower_pos_A, upper_pos_A, lower_pos_B, upper_pos_B);
    printf("col %d :: A( %d:%d) B(%d:%d)\n", a_column_pos, lower_limit_A, upper_limit_A, lower_limit_B, upper_limit_B);
    JA1_C[a_column_pos] = c_value_pos;
    for (MKL_INT B_pos = lower_limit_B; B_pos < upper_limit_B; B_pos++){
      for (MKL_INT A_pos = lower_limit_A; A_pos < upper_limit_A; A_pos++, c_value_pos++ ){
        csc_values_C[c_value_pos] = csc_values_B[B_pos] * csc_values_A[A_pos];
        MKL_INT A_padding = IA1_A[A_pos] * number_rows_B;
        IA1_C[c_value_pos] = A_padding;
        printf("\t\t%d ( %d(%f) %d(%f) ) -- %f \n", c_value_pos, A_pos, csc_values_A[A_pos], B_pos , csc_values_B[B_pos], csc_values_C[c_value_pos] );
      }
    }
  }
  /* for (MKL_INT column_krao = 0; column_krao < number_columns; ++column_krao ){
     MKL_INT padding_row = column_krao * number_rows_A;
     MKL_INT padding_c_row = column_krao * number_rows_C;
  // go through all row elements of column X of matrix B
  for ( MKL_INT B_column_pos = 0; B_column_pos < number_rows_B ; ++B_column_pos ){
  // go through all row elements of column X of matrix A

  MKL_INT pos_B = padding_row + B_column_pos;
  for (MKL_INT A_column_pos = 0; A_column_pos < number_rows_A; ++A_column_pos, ++c_value_pos ){
  MKL_INT pos_A = padding_row + A_column_pos;
  MKL_INT pos_C = padding_c_row + ( B_column_pos * number_rows_B) + A_column_pos;
  csc_values_C[pos_C] = csc_values_B[pos_B] * csc_values_A[pos_A];
  printf("%d %d ( %d(%f) %d(%f) ) -- %f \n", c_value_pos, pos_C, pos_A, csc_values_A[pos_A], pos_B,csc_values_B[pos_B], csc_values_C[pos_C] );
  }
  }
  }*/

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
  return 0;
}