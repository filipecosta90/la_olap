#include <stdio.h>
#include <glib.h>
#include <stdlib.h>
#include <string.h>
#include "mkl.h"

//Cache-Lines size is (typically) 64 bytes
#define MEM_LINE_SIZE 64
#define ARRAY_SIZE MEM_LINE_SIZE / sizeof (int)
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

void check_errors(int stat ){
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
  int* aux_coo_rows;
  aux_coo_rows = (int*) malloc (current_values_size * sizeof(int));

  FILE* stream = fopen(argv[1], "r");
  char line[1024];
  int column = atoi (argv[2]);
  MKL_INT number_rows = - 1;
  MKL_INT number_columns = -1 ;
  MKL_INT element_number = 1;
  
    for( element_number = 0 ; (fgets(line, MAX_REG_SIZE, stream) ) ; ++element_number )
  {
    char *field = (char*) malloc( MAX_FIELD_SIZE * sizeof(char) );
    char* tmp_field = strdup(line);
    field = getfield(tmp_field, column, field);
    unsigned int quark_field;
    quark_field = g_quark_from_string (field);

      //printf("%u %d\n", quark_field, number_rowIntel MKL ERROR: Parameter 1 was incorrect on entry to MKL_SCSRCOO.s);
    /* if arrays are full double its size */
    if ( element_number > current_values_size ){

      current_values_size *= GROWTH_FACTOR;
      aux_coo_rows = (int*) realloc(aux_coo_rows, (current_values_size) * GROWTH_FACTOR * sizeof(MKL_INT) );

    }

    /* normal coo property */
    aux_coo_rows[element_number]= quark_field;

    if ( (short) quark_field > number_rows ) {
      number_rows = ( short ) quark_field;
    }
      
    free(tmp_field);
  }
    MKL_INT NNZ = element_number - 1 ;
    number_columns = element_number;

    //define COO sparse-matrix M
    float* coo_values;
    MKL_INT* coo_rows;
    MKL_INT* coo_columns;
    
    coo_values = mkl_malloc ((element_number * sizeof(float)), MEM_LINE_SIZE );
    coo_rows =  mkl_malloc ((element_number * sizeof(MKL_INT)), MEM_LINE_SIZE );
    coo_columns =  mkl_malloc ((element_number * sizeof(MKL_INT)), MEM_LINE_SIZE );

    for (MKL_INT pos = 0; pos < element_number; ++pos) {
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
  
  float* csr_values = NULL;
    long long job[6];
    job[0]=2;
    job[1]=0;
    job[2]=0;
    job[3]=0;
  job[4]=NNZ;
    job[5]=0;

    MKL_INT* AJ;
    MKL_INT* AI;
    csr_values = mkl_malloc ((element_number * sizeof(float)), MEM_LINE_SIZE );
    AJ = mkl_malloc ((element_number * sizeof(float)), MEM_LINE_SIZE );
    AI = mkl_malloc ((element_number+1 * sizeof(float)), MEM_LINE_SIZE );
    
    printf("going to convert matrix N x M = %d x %d\n", number_rows, number_columns);
    int info;
    mkl_scsrcoo (&job[0], &number_columns, csr_values, AJ, AI, &NNZ, coo_values, coo_rows, coo_columns, &info);
    
    check_errors(info);

    return 0;
}
