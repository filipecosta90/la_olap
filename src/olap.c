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
#define ARRAY_SIZE MEM_LINE_SIZE / sizeof (int)
#define GROWTH_FACTOR 2
#define MAX_FIELD_SIZE 128
#define MAX_REG_SIZE 1024

#define MBLK  2

char* getfield(char* line, int num, char* return_string ){
  return_string = strtok(line, "|\n");
  int pos = 1;
  for ( ; pos < num; pos++ ){
    return_string = strtok(NULL, "|\n");
  }
  return return_string;
}

void check_errors(MKL_INT stat ){
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
  printf("hello\n");
  printf("going to convert matrix N x M = %d x %d\nNNZ: %d JOB4: %d\n", number_rows, number_columns, NNZ, job[4]);
  MKL_INT info=-1;
  mkl_scsrcoo (job, &number_rows, csr_values, JA, IA, &NNZ, coo_values, coo_rows, coo_columns, &info);
  check_errors(info);

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
  //   CONVERT FROM CSR TO BSR
  //
  ////////////////////////////////

  job[0] = 0;  //If job[0]=0, the matrix in the CSR format is converted to the BSR format;
  job[1] = 0;  //If job[1]=0, zero-based indexing for the matrix in CSR format is used;
  job[2] = 0;  //If job[2]=0, zero-based indexing for the matrix in the BSR format is used;
  job[3] = 0;  //
  job[4] = 0;  //
  job[5] = 1;  //If job[5]>0, all output arrays absr, jab, and iab are filled in for the
  // output storage.

  float* bsr_values = NULL;
  MKL_INT* JAB;
  MKL_INT* IAB;
  MKL_INT mblk = MBLK;
  MKL_INT ldAbsr;
  ldAbsr=mblk*mblk;
  MKL_INT NROWS =  ( number_rows / mblk ) +1 ;
  bsr_values = (float*) mkl_malloc ((NNZ * sizeof(float)), MEM_LINE_SIZE );
  JAB = (MKL_INT*) mkl_malloc( ((MBLK+1) * sizeof(MKL_INT)), MEM_LINE_SIZE );
  IAB = (MKL_INT*) mkl_malloc (((NROWS+1) * sizeof(MKL_INT)), MEM_LINE_SIZE );

  MKL_INT info_bsr=-1;

  for ( MKL_INT i=0; i < NNZ; i++){
    csr_values[i]=0.0;
  }

  printf("going to convert to BSR:\n");
  mkl_scsrbsr ( job ,&number_rows, &mblk, &ldAbsr, csr_values, JA, IA, bsr_values, JAB, IAB, &info_bsr );
  check_errors(info_bsr);

  for ( int i=0; i < NNZ; i++){
    printf("%f, ", bsr_values[i]);
  }
  printf("\n");
  for ( int i=0; i < mblk; i++){ 
    printf("%d, ", JAB[i]);
  }
  printf("\n");
  for ( int i=0; i < NROWS+1; i++){
    printf("%d, ", IAB[i]);
  }
  printf("\n");
  return 0;
}
