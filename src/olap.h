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

void tbl_read( char* table_name, MKL_INT tbl_column, MKL_INT* nnz, MKL_INT* rows, MKL_INT* columns , float* A_csr_values, MKL_INT* A_JA, MKL_INT* A_IA){
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

  for( element_number = 0 ; (fgets(line, MAX_REG_SIZE, stream) ) ; ++element_number )
  {
    char *field = (char*) malloc( MAX_FIELD_SIZE * sizeof(char) );
    char* tmp_field = strdup(line);
    field = getfield(tmp_field, tbl_column, field);
    MKL_INT quark_field;

    quark_field = g_quark_from_string (field);

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

  A_csr_values = (float*) mkl_malloc ((element_number * sizeof(float)), MEM_LINE_SIZE );
  A_JA = (MKL_INT*) mkl_malloc (( element_number * sizeof(MKL_INT)), MEM_LINE_SIZE );
  A_IA = (MKL_INT*) mkl_malloc ((number_rows+1 * sizeof(MKL_INT)), MEM_LINE_SIZE );

  sparse_status_t status_coo_csr;
  mkl_scsrcoo (job, &number_rows, A_csr_values, A_JA, A_IA, &NNZ, coo_values, coo_rows, coo_columns, &status_coo_csr);
  mkl_free(coo_values);
  mkl_free(coo_rows);
  mkl_free(coo_columns);
}

