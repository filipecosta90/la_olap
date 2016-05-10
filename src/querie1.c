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

  //define CSR sparse-matrix M

  float* returnFlag_csr_values = NULL;
  MKL_INT* returnFlag_JA;
  MKL_INT* returnFlag_IA;
  MKL_INT returnFlag_rows;
  MKL_INT returnFlag_columns;
  MKL_INT returnFlag_nnz;

  float* lineStatus_csr_values = NULL;
  MKL_INT* lineStatus_JA;
  MKL_INT* lineStatus_IA;
  MKL_INT lineStatus_rows;
  MKL_INT lineStatus_columns;
  MKL_INT lineStatus_nnz;

  float* quantity_csr_values = NULL;
  MKL_INT* quantity_JA;
  MKL_INT* quantity_IA;
  MKL_INT quantity_rows;
  MKL_INT quantity_columns;
  MKL_INT quantity_nnz;

  //read return flag
  tbl_read( "__tbl/lineitem.tbl" , 9, &returnFlag_nnz, &returnFlag_rows, &returnFlag_columns , &returnFlag_csr_values, &returnFlag_JA, &returnFlag_IA);

  //read line status
  tbl_read( "__tbl/lineitem.tbl" , 10, &lineStatus_nnz, &lineStatus_rows, &lineStatus_columns , &lineStatus_csr_values, &lineStatus_JA, &lineStatus_IA);

  //read quantity
  tbl_read( "__tbl/lineitem.tbl" , 4, &quantity_nnz, &quantity_rows, &quantity_columns , &quantity_csr_values, &quantity_JA, &quantity_IA);


  // compute C = A krao B
  /*  csr_kron( 
      A_csr_values, A_JA, A_IA, A_nnz, A_rows, A_columns, 
      B_csr_values, B_JA, B_IA, B_nnz, A_rows, B_columns, 
      &C_csr_values, &C_JA, &C_IA, &C_nnz, &C_rows, &C_columns 
      );
      */  
  return 0;
}
