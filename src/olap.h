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

//starts at position 1 
char* getfield(char* line, int num, char* return_string ){
  return_string = strtok(line, "|\n");
  int pos = 1;
  for ( ; pos < num; pos++ ){
    return_string = strtok(NULL, "|\n");
  }
  return return_string;
}

void print_csc(float* csc_values, MKL_INT* JA1, MKL_INT* IA1, MKL_INT NNZ, MKL_INT number_rows, MKL_INT number_columns ){
  printf("N NONZ: %d\t", NNZ);
  printf("N ROWS: %d\t", number_rows);
  printf("N COLS: %d\n", number_columns);
  printf("CSC VALUES(%d):\t", sizeof(csc_values));
  for (MKL_INT pos = 0; pos < NNZ; pos++){
    printf("%f, ", csc_values[pos]);
  }
  printf("\nJA1:\t");

  for (int pos = 0; pos < NNZ; pos++){
    printf("%d, ", JA1[pos]);
  }

  printf("\nIA1:\t");
  for (int pos = 0; pos <= number_columns; pos++){
    printf("%d, ", IA1[pos]);
  }
  printf("\n");
}

void print_csr(float* csr_values, MKL_INT* JA, MKL_INT* IA, MKL_INT NNZ, MKL_INT number_rows, MKL_INT number_columns ){
  printf("N NONZ: %d\t", NNZ);
  printf("N ROWS: %d\t", number_rows);
  printf("N COLS: %d\n", number_columns);
  printf("CSR VALUES(%d):\t", sizeof(csr_values));
  for (MKL_INT pos = 0; pos < NNZ; pos++){
    printf("%f, ", csr_values[pos]);
  }
  printf("\nJA:\t");

  for (int pos = 0; pos < NNZ; pos++){
    printf("%d, ", JA[pos]);
  }

  printf("\nIA:\t");
  for (int pos = 0; pos <= number_rows; pos++){
    printf("%d, ", IA[pos]);
  }
  printf("\n");
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

void tbl_read( char* table_name, MKL_INT tbl_column, MKL_INT* nnz, MKL_INT* rows, MKL_INT* columns , float** A_csr_values, MKL_INT** A_JA, MKL_INT** A_IA){
  MKL_INT current_values_size = ARRAY_SIZE;
  //define COO sparse-matrix M
  MKL_INT* aux_coo_rows;
  aux_coo_rows = (MKL_INT*) malloc (current_values_size * sizeof(MKL_INT));

  FILE* stream = fopen(table_name, "r");
  char line[1024];
  MKL_INT number_rows = - 1;
  MKL_INT number_columns = -1 ;
  MKL_INT element_number = 1;
  MKL_INT job[8];
  MKL_INT padding_quark = 0;
  for( element_number = 0 ; (fgets(line, MAX_REG_SIZE, stream) ) ; ++element_number )
  {
    char *field = (char*) malloc( MAX_FIELD_SIZE * sizeof(char) );
    char* tmp_field = strdup(line);
    field = getfield(tmp_field, tbl_column, field);
    MKL_INT quark_field;

    quark_field = g_quark_from_string (field);
    if (quark_field > 1 && element_number == 0 ){
      padding_quark = quark_field - 1;
    }
    quark_field -= padding_quark;
    /* if arrays are full double its size */
    if ( element_number >= current_values_size ){
      current_values_size *= GROWTH_FACTOR;
      aux_coo_rows = (MKL_INT*) realloc(aux_coo_rows, (current_values_size) * GROWTH_FACTOR * sizeof(MKL_INT) );
    }

    /* normal coo property */
    aux_coo_rows[element_number]=  quark_field - 1 ;

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
  }
  free(aux_coo_rows);

  /////////////////////////////////
  //
  //   CONVERT FROM COO TO CSR
  //
  ////////////////////////////////

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

  *A_csr_values = (float*) mkl_malloc ((NNZ * sizeof(float)), MEM_LINE_SIZE );
  *A_JA = (MKL_INT*) mkl_malloc (( NNZ * sizeof(MKL_INT)), MEM_LINE_SIZE );
  *A_IA = (MKL_INT*) mkl_malloc ((number_rows+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );

  sparse_status_t status_coo_csr;
  mkl_scsrcoo (job, &number_rows, *A_csr_values, *A_JA, *A_IA, &NNZ, coo_values, coo_rows, coo_columns, &status_coo_csr);
  check_errors(status_coo_csr);
  mkl_free(coo_values);
  mkl_free(coo_rows);
  mkl_free(coo_columns);
  *rows = number_rows;
  *columns = number_columns;
  *nnz = NNZ;
}

/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
//
//   COMPUTE HADAMARD
//
/////////////////////////////////
void csr_hadamard(
    float* A_csr_values, MKL_INT* A_JA, MKL_INT* A_IA, MKL_INT A_NNZ, MKL_INT number_rows, 
    float* B_csr_values, MKL_INT* B_JA, MKL_INT* B_IA , MKL_INT B_NNZ,
    float** C_csr_values, MKL_INT** C_JA, MKL_INT** C_IA, MKL_INT *C_NNZ
    ){
  MKL_INT NNZ; 
  if (A_NNZ > B_NNZ ){
    NNZ = A_NNZ;
  }
  else {
    NNZ = B_NNZ;
  }
  *C_csr_values = (float*) mkl_malloc ((NNZ * sizeof(float)), MEM_LINE_SIZE );
  *C_JA = (MKL_INT*) mkl_malloc (( NNZ * sizeof(MKL_INT)), MEM_LINE_SIZE );
  *C_IA = (MKL_INT*) mkl_malloc ((number_rows+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );

  MKL_INT c_pos = 0;
  MKL_INT at_row = 0;
  for ( ; at_row < number_rows; ++at_row){
    // insert start of line int C_IA
    (*C_IA)[at_row] = c_pos;
    //pivot positions
    MKL_INT column_A_pivot = A_IA[at_row];
    MKL_INT column_B_pivot = B_IA[at_row];
    //limit positions
    MKL_INT column_A_limit = A_IA[at_row+1];
    MKL_INT column_B_limit = B_IA[at_row+1];

    MKL_INT A_line_sizeof = column_A_limit - column_A_pivot;
    MKL_INT B_line_sizeof = column_B_limit - column_B_pivot;

    if (A_line_sizeof < B_line_sizeof){
      for ( ; column_A_pivot < column_A_limit  ; ++column_A_pivot ){
        for ( ; A_JA[column_A_pivot] < B_JA[column_B_pivot] && column_B_pivot < column_B_limit ; ++column_B_pivot){
        }
        if ( A_JA[column_A_pivot] == B_JA[column_B_pivot] ){
          //insert into C
          (*C_csr_values)[c_pos] = A_csr_values[c_pos] * B_csr_values[c_pos];
          (*C_JA)[c_pos]=column_A_pivot;
          ++c_pos;
          ++column_B_pivot;
        }

      }
    }
    else {
      for ( ; column_B_pivot < column_B_limit  ; ++column_B_pivot){
        for ( ; B_JA[column_A_pivot] <  A_JA[column_B_pivot] && column_A_pivot < column_A_limit; ++column_A_pivot ){
        }
        if ( A_JA[column_A_pivot] == B_JA[column_B_pivot] ){
          //insert into C
          (*C_csr_values)[c_pos] = A_csr_values[c_pos] * B_csr_values[c_pos];
          (*C_JA)[c_pos]=column_B_pivot;
          ++column_A_pivot;
          ++c_pos;
        }
      }
    }
  }
  //insert the final C_JA position 
  (*C_IA)[at_row]=c_pos;
  *C_NNZ = c_pos;
}

/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
//
//   COMPUTE KHATRI-RAO
//
/////////////////////////////////
void csr_krao( 
    float* A_csr_values, MKL_INT* A_JA, MKL_INT* A_IA, MKL_INT A_NNZ, MKL_INT A_number_rows, MKL_INT A_number_columns,
    float* B_csr_values, MKL_INT* B_JA, MKL_INT* B_IA , MKL_INT B_NNZ, MKL_INT B_number_rows, MKL_INT B_number_columns,
    float** C_csr_values, MKL_INT** C_JA, MKL_INT** C_IA, MKL_INT* C_NNZ, MKL_INT* C_number_rows, MKL_INT* C_number_columns  
    ){

  MKL_INT job[8];
  /////////////////////////////////////
  //
  //   CONVERT A and B from CSR to CSC
  //
  /////////////////////////////////////
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
  sparse_status_t status_convert_csc;
  float* A_csc_values = NULL;
  MKL_INT* A_JA1;
  MKL_INT* A_IA1;

  A_csc_values = (float*) mkl_malloc (( A_NNZ * sizeof(float) ), MEM_LINE_SIZE );
  A_JA1 = (MKL_INT*) mkl_malloc (( A_NNZ * sizeof(MKL_INT) ), MEM_LINE_SIZE );
  A_IA1 = (MKL_INT*) mkl_malloc ((A_number_columns+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );
  mkl_scsrcsc(job, &A_NNZ, A_csr_values, A_JA, A_IA, A_csc_values, A_JA1, A_IA1, &status_convert_csc);
  printf("A csr -> csc ok?\n");
  check_errors(status_convert_csc);

  print_csc(A_csc_values, A_JA1, A_IA1, A_NNZ, A_number_rows, A_number_columns );
  float* B_csc_values = NULL;
  MKL_INT* B_JA1;
  MKL_INT* B_IA1;

  B_csc_values = (float*) mkl_malloc (( B_NNZ * sizeof(float)), MEM_LINE_SIZE );
  B_JA1 = (MKL_INT*) mkl_malloc (( B_NNZ * sizeof(MKL_INT)), MEM_LINE_SIZE );
  B_IA1 = (MKL_INT*) mkl_malloc (( B_number_columns+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );
  mkl_scsrcsc(job, &B_NNZ, B_csr_values, B_JA, B_IA, B_csc_values, B_JA1, B_IA1, &status_convert_csc);
  printf("A csr -> csc ok?\n");
  check_errors(status_convert_csc);
  print_csc(B_csc_values, B_JA1, B_IA1, B_NNZ, B_number_rows, B_number_columns );
  /////////////////////////////////
  //
  //   COMPUTE KRAO
  //
  ////////////////////////////////
  float* C_csc_values = NULL;
  MKL_INT* C_JA1;
  MKL_INT* C_IA1;

  C_csc_values = (float*) mkl_malloc (( A_NNZ * sizeof(float)), MEM_LINE_SIZE );
  C_JA1 = (MKL_INT*) mkl_malloc (( A_NNZ  * sizeof(MKL_INT)), MEM_LINE_SIZE );
  C_IA1 = (MKL_INT*) mkl_malloc (( A_number_columns+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );

  MKL_INT at_column = 0;
  MKL_INT max_row = 0;
  MKL_INT row_pos = 0;

  for ( ; at_column < A_number_columns; ++at_column){
    // insert start of column int C_IA1
    C_IA1[at_column] = A_IA1[at_column];
    C_csc_values[at_column] = A_csc_values[at_column] * B_csc_values[at_column];
    row_pos  = ( A_JA1[at_column] * B_number_rows ) + B_JA1[at_column];
    C_JA1[at_column] = row_pos;
    if (max_row < row_pos){
      max_row = row_pos;
    }
  }
  ++max_row;
  C_IA1[at_column] = A_NNZ;
  print_csc(C_csc_values, C_JA1, C_IA1, A_NNZ, max_row, A_number_columns );
  /////////////////////////////////
  //
  //   CONVERT C from CSC to CSR
  //
  ////////////////////////////////

  // If job[0]=1, the matrix in the CSC format is converted to the CSR format.
  job[0] = 1;

  // job[1]
  // If job[1]=0, zero-based indexing for the matrix in CSR format is used;
  // if job[1]=1, one-based indexing for the matrix in CSR format is used.
  job[1] = 0;

  // job[2]
  // If job[2]=0, zero-based indexing for the matrix in the CSC format is used;
  // if job[2]=1, one-based indexing for the matrix in the CSC format is used.
  job[2] = 0;

  // job[5] - job indicator.
  // If job[5]=0, only arrays ja1, ia1 are filled in for the output storage.
  // If job[5]≠0, all output arrays acsc, ja1, and ia1 are filled in for the output storage.
  job[5] = 1;
  sparse_status_t status_convert_csr;

  *C_csr_values = (float*) mkl_malloc (( A_NNZ * sizeof(float)), MEM_LINE_SIZE );
  *C_JA = (MKL_INT*) mkl_malloc (( A_NNZ * sizeof(MKL_INT)), MEM_LINE_SIZE );
  *C_IA = (MKL_INT*) mkl_malloc (( max_row + 1 * sizeof(MKL_INT)), MEM_LINE_SIZE );
  mkl_scsrcsc(job, &A_NNZ, *C_csr_values, *C_JA, *C_IA, C_csc_values, C_JA1, C_IA1, &status_convert_csr);
  check_errors(status_convert_csr);
  printf("C csc -> csr ok?\n");
  *C_number_rows = max_row ;
  *C_number_columns = A_number_columns;
  *C_NNZ = A_NNZ;
}


/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
/////////////////////////////////
//
//   COMPUTE KRONECKER PRODUCT
//
/////////////////////////////////

void csr_kron(
    float* A_csr_values, MKL_INT* A_JA, MKL_INT* A_IA, MKL_INT A_NNZ, MKL_INT A_number_rows, MKL_INT A_number_columns,
    float* B_csr_values, MKL_INT* B_JA, MKL_INT* B_IA , MKL_INT B_NNZ, MKL_INT B_number_rows, MKL_INT B_number_columns,
    float** C_csr_values, MKL_INT** C_JA, MKL_INT** C_IA, MKL_INT* C_NNZ, MKL_INT* C_number_rows, MKL_INT* C_number_columns
    ){

  MKL_INT job[8];
  /////////////////////////////////////
  //
  //   CONVERT A and B from CSR to CSC
  //
  /////////////////////////////////////
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
  sparse_status_t status_convert_csc;
  float* A_csc_values = NULL;
  MKL_INT* A_JA1;
  MKL_INT* A_IA1;

  A_csc_values = (float*) mkl_malloc (( A_NNZ * sizeof(float) ), MEM_LINE_SIZE );
  A_JA1 = (MKL_INT*) mkl_malloc (( A_NNZ * sizeof(MKL_INT) ), MEM_LINE_SIZE );
  A_IA1 = (MKL_INT*) mkl_malloc ((A_number_columns+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );
  mkl_scsrcsc(job, &A_NNZ, A_csr_values, A_JA, A_IA, A_csc_values, A_JA1, A_IA1, &status_convert_csc);
  printf("A csr -> csc ok?\n");
  check_errors(status_convert_csc);

  print_csc(A_csc_values, A_JA1, A_IA1, A_NNZ, A_number_rows, A_number_columns );
  float* B_csc_values = NULL;
  MKL_INT* B_JA1;
  MKL_INT* B_IA1;

  B_csc_values = (float*) mkl_malloc (( B_NNZ * sizeof(float)), MEM_LINE_SIZE );
  B_JA1 = (MKL_INT*) mkl_malloc (( B_NNZ * sizeof(MKL_INT)), MEM_LINE_SIZE );
  B_IA1 = (MKL_INT*) mkl_malloc (( B_number_columns+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );
  mkl_scsrcsc(job, &B_NNZ, B_csr_values, B_JA, B_IA, B_csc_values, B_JA1, B_IA1, &status_convert_csc);
  printf("B csr -> csc ok?\n");
  check_errors(status_convert_csc);
  print_csc(B_csc_values, B_JA1, B_IA1, B_NNZ, B_number_rows, B_number_columns );
  /////////////////////////////////
  //
  //   COMPUTE KRON
  //
  ////////////////////////////////
  float* C_csc_values = NULL;
  MKL_INT* C_JA1;
  MKL_INT* C_IA1;
  MKL_INT C_nnz = A_NNZ * B_NNZ;
  MKL_INT C_ncols = A_number_columns * B_number_columns;
  C_csc_values = (float*) mkl_malloc (( C_nnz * sizeof(float)), MEM_LINE_SIZE );
  C_JA1 = (MKL_INT*) mkl_malloc (( C_nnz * sizeof(MKL_INT)), MEM_LINE_SIZE );
  C_IA1 = (MKL_INT*) mkl_malloc (( C_ncols +1 * sizeof(MKL_INT)), MEM_LINE_SIZE );

  MKL_INT at_column = 0;
  MKL_INT at_column_A = 0;
  MKL_INT at_column_B = 0;
  MKL_INT max_row = 0;
  MKL_INT row_pos = 0;


  for ( ; at_column_A < A_number_columns; ++at_column_A ){
    at_column_B = 0;
    for ( ; at_column_B < B_number_columns; ++at_column_B ){
      at_column = A_IA1[at_column_A] * B_number_columns + B_IA1[at_column_B];
      C_IA1[at_column] = at_column;
      C_csc_values[at_column] = A_csc_values[at_column_A] * B_csc_values[at_column_B];
      row_pos  = ( A_JA1[at_column_A] * B_number_rows ) + B_JA1[at_column_B];
      C_JA1[at_column] = row_pos;
      if (max_row < row_pos){
        max_row = row_pos;
      }
    }
  }
  at_column++;
  ++max_row;
  C_IA1[at_column] = C_nnz;
  print_csc(C_csc_values, C_JA1, C_IA1, C_nnz, max_row,   C_ncols );
  printf("printed\n");
  /////////////////////////////////
  //
  //   CONVERT C from CSC to CSR
  //
  ////////////////////////////////

  // If job[0]=1, the matrix in the CSC format is converted to the CSR format.
  job[0] = 1;

  // job[1]
  // If job[1]=0, zero-based indexing for the matrix in CSR format is used;
  // if job[1]=1, one-based indexing for the matrix in CSR format is used.
  job[1] = 0;

  // job[2]
  // If job[2]=0, zero-based indexing for the matrix in the CSC format is used;
  // if job[2]=1, one-based indexing for the matrix in the CSC format is used.
  job[2] = 0;

  // job[5] - job indicator.
  // If job[5]=0, only arrays ja1, ia1 are filled in for the output storage.
  // If job[5]≠0, all output arrays acsc, ja1, and ia1 are filled in for the output storage.
  job[5] = 1;
  sparse_status_t status_convert_csr;

  *C_csr_values = (float*) mkl_malloc (( C_nnz * sizeof(float)), MEM_LINE_SIZE );
  *C_JA = (MKL_INT*) mkl_malloc (( C_nnz * sizeof(MKL_INT)), MEM_LINE_SIZE );
  *C_IA = (MKL_INT*) mkl_malloc (( max_row + 1 * sizeof(MKL_INT)), MEM_LINE_SIZE );
  printf("going to convert back to CSR\n");
  mkl_scsrcsc(job, &C_nnz, *C_csr_values, *C_JA, *C_IA, C_csc_values, C_JA1, C_IA1, &status_convert_csr);
  check_errors(status_convert_csr);
  printf("C csc -> csr ok?\n");
  *C_number_rows = max_row ;
  *C_number_columns = C_ncols;
  *C_NNZ = C_nnz;
}


