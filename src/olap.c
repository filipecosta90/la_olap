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
#define MAX_FIELD_SIZE 128
#define MAX_REG_SIZE 1024

char* getfield(char* line, int num, char* return_string ){
  return_string = strtok(line, "|");
  int pos = 1;
  for ( ; pos <= num; pos++ ){
    return_string = strtok(NULL, "|\n");
  }
  return return_string;
}

int main( int argc, char* argv[]){

  int current_values_size = ARRAY_SIZE;

  //define COO sparse-matrix M
  float* coo_values;
  int* coo_rows;
  int* coo_columns;
  coo_values = (float*) malloc (current_values_size * sizeof(float));
  coo_rows = (int*) malloc (current_values_size * sizeof(int));
  coo_columns = (int*) malloc (current_values_size * sizeof(int));

  FILE* stream = fopen(argv[1], "r");
  char line[1024];
  int column = atoi (argv[2]);
  int max_row = -1;
  int max_column = -1;
  int element_number = 0;
  for( element_number = 0 ; (fgets(line, MAX_REG_SIZE, stream) ) ; ++element_number )
  {
    char *field = (char*) malloc( MAX_FIELD_SIZE * sizeof(char) );
    char* tmp_field = strdup(line);
    field = getfield(tmp_field, column, field);
    GQuark* quark_field;
    quark_field = g_quark_from_string (field);

    /* if arrays are full double its size */
    if ( element_number > current_values_size ){

      current_values_size *= GROWTH_FACTOR;
      coo_values = (float*) realloc(coo_values, (current_values_size+1) * GROWTH_FACTOR * sizeof(float) );
      coo_rows = (int*) realloc(coo_rows, (current_values_size+1) * GROWTH_FACTOR * sizeof(int) );
      coo_columns = (int*) realloc(coo_columns, (current_values_size+1) * GROWTH_FACTOR * sizeof(int) );

    }

    /* particular property */
    coo_values[element_number]=1.0;
    coo_rows[element_number]=column;

    /* normal csr property */
    coo_columns[element_number]=quark_field;

    if ((int) quark_field > max_row) {
      max_row =  quark_field;
    }

    free(tmp_field);
  }

  coo_values = (float*) realloc(coo_values, (element_number+1) * GROWTH_FACTOR * sizeof(float) );
  coo_rows = (int*) realloc(coo_rows, (element_number+1) * GROWTH_FACTOR * sizeof(int) );
  coo_columns = (int*) realloc(coo_columns, (element_number+1) * GROWTH_FACTOR * sizeof(int) );

  sparse_matrix_t table_projection_coo;
  int stat = -1 ;
  stat = mkl_sparse_s_create_coo(&table_projection_coo, SPARSE_INDEX_BASE_ZERO ,  max_row, element_number,  element_number, coo_rows, coo_columns, coo_values );
  if ( stat == SPARSE_STATUS_SUCCESS){
    printf( "SPARSE_STATUS_SUCCESS\n");
  }

  /////////////////////////////////
  //
  //   CONVERT FROM COO TO CSR
  //
  ////////////////////////////////

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
  job[3]=0;

  //job[4]=nzmax - maximum number of the non-zero elements allowed if job[0]=0.
  job[4]=element_number;

  // For conversion to the CSR format:
  // If job[5]=0, all arrays acsr, ja, ia are filled in for the output storage.
  // If job[5]=1, only array ia is filled in for the output storage.
  // If job[5]=2, then it is assumed that the routine already has been called with the job[5]=1,
  // and the user allocated the required space for storing the output arrays acsr and ja.
  job[5]=1;

  //flag containing the convert result
  int convert_stat;
  float* csr_values;
  int* csr_ia;
  int* csr_ja;
  csr_values = (float*) malloc (element_number+2 * sizeof(float));
  csr_ja = (int*) malloc (element_number+2 * sizeof(int));
  csr_ia = (int*) malloc (element_number+2 * sizeof(int));
  int n = element_numberyt;

  mkl_scsrcoo (&job , &n , csr_values , csr_ja , csr_ia , &element_number , coo_values , coo_rows , coo_columns , &convert_stat );

  if ( convert_stat == 0){
    printf( "The convert execution was successful.\n");
  }


  /*

  //If job[0]=0, the matrix in the CSR format is converted to the BSR format;
  //if job[0]=1, the matrix in the BSR format is converted to the CSR format.
  job[0]=0;

  //If job[1]=0, zero-based indexing for the matrix in CSR format is used;
  //if job[1]=1, one-based indexing for the matrix in CSR format is used.
  job[1]=0;

  //If job[2]=0, zero-based indexing for the matrix in the BSR format is used;
  //if job[2]=1, one-based indexing for the matrix in the BSR format is used.
  job[2]=0;

  //job[3] is only used for conversion to CSR format. By default, the converter 
  //saves the blocks without checking whether an element is zero or not. 
  //If job[3]=1, then the converter only saves non-zero elements in blocks.
  job[3]=0;  
  job[4]=0;  
  //job[5] - job indicator.
  //For conversion to the BSR format:
  //If job[5]=0, only arrays jab, iab are generated for the output storage.
  //If job[5]>0, all output arrays absr, jab, and iab are filled in for the output storage.
  //If job[5]=-1, iab[m] - iab[0] returns the number of non-zero blocks.
  job[5]=1;

  // m 
  // Actual row dimension of the matrix A for convert to the BSR format; 
  // block row dimension of the matrix A for convert to the CSR format.
  // in our case m = max_row;

  // mblk
  // Size of the block in the matrix A.

  // ldabsr
  // Leading dimension of the array absr as declared in the calling program. 
  // ldabsr must be greater than or equal to mblk*mblk.

  // acsr
  // (input/output)
  // Array containing non-zero elements of the matrix A. 
  // Its length is equal to the number of non-zero elements in the matrix A. 
  // Refer to values array description in Sparse Matrix Storage Formats for more details.

  // ja
  // (input/output). Array containing the column indices for each non-zero element of the matrix A.
  // Its length is equal to the length of the array acsr. 
  // Refer to columns array description in Sparse Matrix Storage Formats for more details.

  // ia
  // (input/output). Array of length m + 1, containing indices of elements in the array acsr, 
  // such that ia[I]] - iab[0] is the index in the array acsr of the first non-zero element from the row I. 
  // The value of ia[m]] - iab[0] is equal to the number of non-zeros. 
  // Refer to rowIndex array description in Sparse Matrix Storage Formats for more details.

  // absr
  // (input/output)
  // Array containing elements of non-zero blocks of the matrix A. 
  // Its length is equal to the number of non-zero blocks in the matrix A multiplied by mblk*mblk. 
  // Refer to values array description in BSR Format for more details.

  // jab
  // (input/output). Array containing the column indices for each non-zero block of the matrix A.
  // Its length is equal to the number of non-zero blocks of the matrix A.
  // Refer to columns array description in BSR Format for more details.

  // iab
  // (input/output). Array of length (m + 1), containing indices of blocks in the array absr, 
  // such that iab[i] - iab[0] is the index in the array absr of the first non-zero element from the i-th row . 
  // The value of iab[m] is equal to the number of non-zero blocks. 
  // Refer to rowIndex array description in BSR Format for more details.

  ///////////////////
  //Output Parameters
  ///////////////////

  //info
  //Integer info indicator only for converting the matrix A from the CSR format.
  //If info=0, the execution is successful.
  //If info=1, it means that mblk is equal to 0.
  //If info=2, it means that ldabsr is less than mblk*mblk and there is no space for all blocks.

  float* bsr_values;
  int* bsr_ia;
  int* bsr_ja;
  bsr_values = (float*) malloc (current_values_size * sizeof(float));
  bsr_ia = (int*) malloc (current_values_size * sizeof(int));
  bsr_ja = (int*) malloc (( element_number +1 ) * sizeof(int));

  MKL_INT mblk = LINE_SIZE;
  MKL_INT ldabsr = mblk*mblk;

  printf("matrix %d >< %d, nonZ %d \n", max_row, element_number, element_number ); 
  printf( "Going to convert from CSR to BSR.\n");
  //  printf("sizeof ia: %d sizeof ja: %d sizeof matrixvalues %d\n", element_number, sizeof(), sizeof(), sizeof() );
  //  mkl_scsrbsr (job , &max_row , &mblk , &ldabsr , matrix_values , pointer_B , matrix_rows , bsr_values , bsr_ia , bsr_ja , &convert_stat );

  if ( convert_stat == 0){
    printf( "The convert execution was successful.\n");
  }
  */
    return 0;
}
