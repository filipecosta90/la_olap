/* ---------------------------------------------------------------------------
 **    Filename: olap.c
 **
 **     License: This file is part of OLAP PROJECT.
 **
 **              OLAP PROJECT is free software: you can redistribute it
 **              and/or modify it under the terms of the GNU General Public
 **              License as published by the Free Software Foundation,
 **              either version 3 of the License, or (at your option)
 **              any later version.
 **
 **              OLAP is distributed in the hope that it will be useful,
 **              but WITHOUT ANY WARRANTY; without even the implied warranty of
 **              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 **              GNU General Public License for more details.
 **
 **              You should have received a copy of the GNU General Public
 **              License along with OLAP.
 **              If not, see <http://www.gnu.org/licenses/>.
 **
 ** Description: The proposed file focus on a typed linear algebra approach,
 **              encoding OLAP functionality solely in terms of Linear Algebra
 **              operations, represented in the CSR and BSR format, and
 **              recurring to Intel MKL and GLIB libraries.
 **              It encodes TPC-H querie-1 into Linear Algebra Operations.
 **
 **
 **     Authors: Filipe Oliveira <a57816@alunos.uminho.pt>
 **          and Sérgio Caldas   <a57779@alunos.uminho.pt>
 **
 ** University of Minho, High Performance Computing Dpt. , April 2016
 ** -------------------------------------------------------------------------*/

#include <stdio.h>
#include <glib.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mkl_spblas.h"
#include "mkl_types.h"
#include "mkl.h"
#include "olap.h"

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
  sparse_matrix_t  quantity_matrix;

  //conversion status from csr arrays into mkl sparse_matrix_t 
  sparse_status_t status_to_csr;

  //read return flag
  tbl_read( "__tbl/lineitem.tbl" , 9, &returnFlag_nnz, &returnFlag_rows, &returnFlag_columns , &returnFlag_csr_values, &returnFlag_JA, &returnFlag_IA);

  //read line status
  tbl_read( "__tbl/lineitem.tbl" , 10, &lineStatus_nnz, &lineStatus_rows, &lineStatus_columns , &lineStatus_csr_values, &lineStatus_JA, &lineStatus_IA);

  //read quantity
  tbl_read_measure( "__tbl/lineitem.tbl" , 5, &quantity_nnz, &quantity_rows, &quantity_columns , &quantity_csr_values, &quantity_JA, &quantity_IA);

  //        convert via sparseBLAS API to Handle containing internal data for 
  //        subsequent Inspector-executor Sparse BLAS operations.
  status_to_csr = mkl_sparse_s_create_csr ( &quantity_matrix , SPARSE_INDEX_BASE_ZERO, 
      quantity_rows, quantity_columns, quantity_IA, quantity_IA+1, quantity_JA, quantity_csr_values );
  printf("quantity convert? :");
  check_errors(status_to_csr);


  float* selection_csr_values = NULL;
  MKL_INT* selection_JA;
  MKL_INT* selection_IA;
  MKL_INT selection_rows;
  MKL_INT selection_columns;
  MKL_INT selection_nnz;
  sparse_matrix_t  selection_matrix;

  //read shipdate gt and lt into selection matrix
  tbl_read_filter_and( "__tbl/lineitem.tbl" , 11, GREATER_EQ , "1998-08-28", LESS_EQ , "1998-12-01", &selection_nnz, &selection_rows, &selection_columns , &selection_csr_values, &selection_JA, &selection_IA);

  //        convert via sparseBLAS API to Handle containing internal data for 
  //        subsequent Inspector-executor Sparse BLAS operations.
  status_to_csr = mkl_sparse_s_create_csr ( &selection_matrix , SPARSE_INDEX_BASE_ZERO, selection_rows, selection_columns, selection_IA, selection_IA+1, selection_JA, selection_csr_values );
  printf("selection convert? :");
  check_errors(status_to_csr);

  float* projection_csr_values = NULL;
  MKL_INT* projection_JA;
  MKL_INT* projection_IA;
  MKL_INT projection_rows;
  MKL_INT projection_columns;
  MKL_INT projection_nnz;
  sparse_matrix_t  projection_matrix;

  float* aggregation_csr_values = NULL;
  MKL_INT* aggregation_JA;
  MKL_INT* aggregation_IA;
  MKL_INT aggregation_rows;
  MKL_INT aggregation_columns;
  MKL_INT aggregation_nnz;

  float* intermediate_csr_values = NULL;
  MKL_INT* intermediate_JA;
  MKL_INT* intermediate_IA;
  MKL_INT intermediate_rows;
  MKL_INT intermediate_columns;
  MKL_INT intermediate_nnz;
  sparse_matrix_t  intermediate_matrix;

  float* final_csr_values = NULL;
  MKL_INT* final_JA;
  MKL_INT* final_IA;
  MKL_INT final_rows;
  MKL_INT final_columns;
  MKL_INT final_nnz;

  // compute projection = returnFlag krao lineStatus
  csr_krao(
      returnFlag_csr_values, returnFlag_JA, returnFlag_IA, 
      returnFlag_nnz, returnFlag_rows, returnFlag_columns,
      lineStatus_csr_values, lineStatus_JA, lineStatus_IA ,
      lineStatus_nnz, lineStatus_rows, lineStatus_columns,
      &projection_csr_values, &projection_JA, &projection_IA, 
      &projection_nnz, &projection_rows, &projection_columns  
      );

  status_to_csr = mkl_sparse_s_create_csr ( &projection_matrix , SPARSE_INDEX_BASE_ZERO, projection_rows, projection_columns, projection_IA, projection_IA+1, projection_JA, projection_csr_values );
  printf("projection convert? :");
  check_errors(status_to_csr);

  // compute aggregation = quantity * bang
  float* bang_vector;
  float* aggregation_vector;
  bang_vector = (float*) mkl_malloc ((quantity_columns * sizeof(float)), MEM_LINE_SIZE );
  aggregation_vector = (float*) mkl_malloc ((quantity_columns * sizeof(float)), MEM_LINE_SIZE );

  sparse_status_t aggregation_result;
  struct matrix_descr descrA;
  descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

  aggregation_result = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, quantity_matrix , descrA, bang_vector, 1.0,  aggregation_vector);
  printf("aggregation result? :");
  check_errors(aggregation_result);

  // compute intermediate_result = projection * selection
  sparse_status_t intermediate_result;

  intermediate_result = mkl_sparse_spmm ( SPARSE_OPERATION_NON_TRANSPOSE , 
      projection_matrix,
      selection_matrix, 
      & intermediate_matrix);
  printf("intermediate result? :");
  check_errors(intermediate_result);

  // compute final_result = intermediate_result * aggregation

  return 0;
}

