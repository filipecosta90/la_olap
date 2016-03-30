#include <stdio.h>
#include <glib.h>
#include <stdlib.h>
#include <string.h>
#include "mkl.h"
#include "mkl_spblas.h"

//Cache-Lines size is (typically) 64 bytes
#define LINE_SIZE 64
#define ARRAY_SIZE LINE_SIZE / sizeof (int)
#define GROWTH_FACTOR 2

char* getfield(char* line, int num, char* return_string ){
  return_string = strtok(line, "|");
  int pos = 1;
  for ( ; pos <= num; pos++ ){
    return_string = strtok(NULL, "|\n");
  }
  return return_string;
}

int main( int argc, char* argv[]){
  //define sparse-matrix M
  float* matrix_values;
  int* matrix_rows;
  int* matrix_columns;
  int current_values_size = ARRAY_SIZE;
  matrix_values = (float*) malloc (current_values_size * sizeof(float));
  matrix_rows = (int*) malloc (current_values_size * sizeof(int));
  matrix_columns = (int*) malloc (current_values_size * sizeof(int));

  FILE* stream = fopen(argv[1], "r");
  char line[1024];

  int key_col = atoi (argv[2]);
  int column = atoi (argv[3]);
  int max_row = -1;
  int max_column = -1;
  int element_number = 0;
  for( element_number = 0 ; (fgets(line, 1024, stream) ) ; ++element_number )
  {
    char *key = (char*) malloc( 128 * sizeof(char) );
    char *field = (char*) malloc( 128 * sizeof(char) );
    char* tmp_key = strdup(line);
    char* tmp_field = strdup(line);
    key = getfield(tmp_key, key_col, key);
    field = getfield(tmp_field, column, field);
    GQuark* quark_key;
    quark_key = g_quark_from_string (key);
    GQuark* quark_field;
    quark_field = g_quark_from_string (field);
    if ( element_number > current_values_size ){
      current_values_size *= GROWTH_FACTOR;
      matrix_values = realloc(matrix_values, current_values_size * GROWTH_FACTOR * sizeof(float) );
      matrix_rows =  realloc(matrix_rows, current_values_size * GROWTH_FACTOR * sizeof(int) );
      matrix_columns =  realloc(matrix_columns, current_values_size * GROWTH_FACTOR * sizeof(int) );
    }

    matrix_values[element_number]=1.0;
    matrix_rows[element_number]=quark_field;
    matrix_columns[element_number]=quark_key;
    if ((int)quark_field > max_row) {
      max_row = quark_field;
    }
    if ((int)quark_key > max_column ){
      max_column = quark_key;
    }
    //printf("%d ,%d, 1\n", (int) quark_field , quark_key );
    free(tmp_key);
    free(tmp_field);

  }

  sparse_matrix_t A;
  int stat;
  stat = mkl_sparse_s_create_coo(&A,SPARSE_INDEX_BASE_ZERO, max_row, max_column, element_number, matrix_rows, matrix_columns, matrix_values );
  if ( stat == SPARSE_STATUS_SUCCESS){
    printf( "SPARSE_STATUS_SUCCESS\n");
  }
  sparse_matrix_t A_csr;

  int job[5];

  // If job[0]=0, the matrix in the CSR format is converted to the coordinate format;
  // if job[0]=1, the matrix in the coordinate format is converted to the CSR format.
  // if job[0]=2, the matrix in the coordinate format is converted to the CSR format, 
  // and the column indices in CSR representation are sorted in the increasing order within each row.
  job[0]=2;

  //If job[1]=0, zero-based indexing for the matrix in CSR format is used;
  //if job[1]=1, one-based indexing for the matrix in CSR format is used.
  job[1]=0;

  //If job[2]=0, zero-based indexing for the matrix in coordinate format is used;
  //if job[2]=1, one-based indexing for the matrix in coordinate format is used.
  job[2]=0;

  //job[4]=nzmax - maximum number of the non-zero elements allowed if job[0]=0.
  job[4]=element_number;

  // For conversion to the CSR format:
  // If job[5]=0, all arrays acsr, ja, ia are filled in for the output storage.
  // If job[5]=1, only array ia is filled in for the output storage.
  // If job[5]=2, then it is assumed that the routine already has been called with the job[5]=1, 
  // and the user allocated the required space for storing the output arrays acsr and ja.
  job[5]=0;

  //if the matrix A was squared what would be its size?
  int matrix_size = max_row > max_column ? max_row : max_column;

  //flag containing the convert result
  int convert_stat;
  float* csr_values;
  int* csr_ia;
  int* csr_ja;
 csr_values = (float*) malloc (current_values_size * sizeof(float));
  csr_ia = (int*) malloc (current_values_size * sizeof(int));
  csr_ja = (int*) malloc (current_values_size * sizeof(int));


  mkl_scsrcoo (&job , &element_number , &csr_values , &csr_ja , &csr_ia , &element_number , &matrix_values , &matrix_rows , &matrix_columns , &convert_stat );
  
   if ( convert_stat == 0){
    printf( "The convert execution was successful.\n");
  }
  
  printf("matrix %d >< %d, nonZ %d \n", max_row, max_column, element_number ); 

  return 0;
}
