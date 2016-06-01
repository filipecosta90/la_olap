/* ---------------------------------------------------------------------------
 **    Filename: querie1_osx.c
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
#include "olap_osx.h"
#include "timer.h"


double global_time_start, global_time_stop, total_time;

void writeResults ( ) {
  total_time = global_time_stop - global_time_start;
  FILE* stream = fopen("timing/timings_osx.dat", "a+");
  fprintf(stream, "1,%f\n", total_time);
  fclose(stream);
}

int main( int argc, char* argv[]){

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

  float* shipdate_gt_csr_values = NULL;
  MKL_INT* shipdate_gt_JA;
  MKL_INT* shipdate_gt_IA;
  MKL_INT shipdate_gt_rows;
  MKL_INT shipdate_gt_columns;
  MKL_INT shipdate_gt_nnz;
  sparse_matrix_t  shipdate_gt_matrix;

  float* shipdate_lt_csr_values = NULL;
  MKL_INT* shipdate_lt_JA;
  MKL_INT* shipdate_lt_IA;
  MKL_INT shipdate_lt_rows;
  MKL_INT shipdate_lt_columns;
  MKL_INT shipdate_lt_nnz;
  sparse_matrix_t  shipdate_lt_matrix;

  float* selection_csr_values = NULL;
  MKL_INT* selection_JA;
  MKL_INT* selection_IA;
  MKL_INT selection_rows;
  MKL_INT selection_columns;
  MKL_INT selection_nnz;
  sparse_matrix_t  selection_matrix;

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
  MKL_INT* final_IA_END;
  MKL_INT final_rows;
  MKL_INT final_columns;
  MKL_INT final_nnz;
  sparse_matrix_t  final_matrix;

  float* bang_vector;
  float* aggregation_vector;

  //conversion status from csr arrays into mkl sparse_matrix_t 
  sparse_status_t status_to_csr;

  //read return flag
  //tbl_read( "__tbl/lineitem.tbl" , 9, &returnFlag_nnz, &returnFlag_rows, &returnFlag_columns, &returnFlag_csr_values, &returnFlag_JA, &returnFlag_IA);
  //convert_and_write_to_csv("__quark_csv/return_flag.mx", returnFlag_csr_values, returnFlag_JA, returnFlag_IA, returnFlag_nnz, returnFlag_rows, returnFlag_columns);
  read_from_mx("__quark_mx_1/return_flag_1.mx", &returnFlag_csr_values, &returnFlag_JA, &returnFlag_IA, &returnFlag_nnz, &returnFlag_rows, &returnFlag_columns);

  //read line status
  //tbl_read( "__tbl/lineitem.tbl" , 10, &lineStatus_nnz, &lineStatus_rows, &lineStatus_columns , &lineStatus_csr_values, &lineStatus_JA, &lineStatus_IA);
  //convert_and_write_to_csv("__quark_csv/line_status.mx", lineStatus_csr_values, lineStatus_JA, lineStatus_IA, lineStatus_nnz, lineStatus_rows, lineStatus_columns);
  read_from_mx("__quark_mx_1/line_status_1.mx", &lineStatus_csr_values, &lineStatus_JA, &lineStatus_IA, &lineStatus_nnz, &lineStatus_rows, &lineStatus_columns);

  //read quantity
  //tbl_read_measure( "__tbl/lineitem.tbl" , 5, &quantity_nnz, &quantity_rows, &quantity_columns , &quantity_csr_values, &quantity_JA, &quantity_IA);
  //convert_and_write_to_csv("__quark_csv/quantity.mx", quantity_csr_values, quantity_JA, quantity_IA, quantity_nnz, quantity_rows, quantity_columns);
  read_from_mx("__quark_mx_1/quantity_1.mx", &quantity_csr_values, &quantity_JA, &quantity_IA, &quantity_nnz, &quantity_rows, &quantity_columns);


  //        convert via sparseBLAS API to Handle containing internal data for 
  //        subsequent Inspector-executor Sparse BLAS operations.
  status_to_csr = mkl_sparse_s_create_csr ( &quantity_matrix , SPARSE_INDEX_BASE_ZERO, 
      quantity_rows, quantity_columns, quantity_IA, quantity_IA+1, quantity_JA, quantity_csr_values );

  //read shipdate gt
  // tbl_read_filter( "__tbl/lineitem.tbl" , 11, GREATER_EQ , "1998-08-28",
  //     &shipdate_gt_nnz, &shipdate_gt_rows, &shipdate_gt_columns , &shipdate_gt_csr_values, &shipdate_gt_JA, &shipdate_gt_IA);
  // convert_and_write_to_csv("__quark_csv/shipdate_gt.mx", shipdate_gt_csr_values, shipdate_gt_JA, shipdate_gt_IA, shipdate_gt_nnz, shipdate_gt_rows, shipdate_gt_columns);
  read_from_mx("__quark_mx_1/shipdate_gt_1.mx", &shipdate_gt_csr_values, &shipdate_gt_JA, &shipdate_gt_IA, &shipdate_gt_nnz, &shipdate_gt_rows, &shipdate_gt_columns);


  //        convert via sparseBLAS API to Handle containing internal data for
  //        subsequent Inspector-executor Sparse BLAS operations.
  status_to_csr = mkl_sparse_s_create_csr ( &shipdate_gt_matrix , SPARSE_INDEX_BASE_ZERO,
      shipdate_gt_rows, shipdate_gt_columns, shipdate_gt_IA, shipdate_gt_IA+1, shipdate_gt_JA, shipdate_gt_csr_values );

  //read shipdate lt
  //  tbl_read_filter( "__tbl/lineitem.tbl" , 11, LESS_EQ , "1998-12-01",
  //    &shipdate_lt_nnz, &shipdate_lt_rows, &shipdate_lt_columns , &shipdate_lt_csr_values, &shipdate_lt_JA, &shipdate_lt_IA);
  //  convert_and_write_to_csv("__quark_csv/shipdate_lt.mx", shipdate_lt_csr_values, shipdate_lt_JA, shipdate_lt_IA, shipdate_lt_nnz, shipdate_lt_rows, shipdate_lt_columns);
  read_from_mx("__quark_mx_1/shipdate_lt_1.mx", &shipdate_lt_csr_values, &shipdate_lt_JA, &shipdate_lt_IA, &shipdate_lt_nnz, &shipdate_lt_rows, &shipdate_lt_columns);

  //        convert via sparseBLAS API to Handle containing internal data for
  //        subsequent Inspector-executor Sparse BLAS operations.
  status_to_csr = mkl_sparse_s_create_csr ( &shipdate_lt_matrix , SPARSE_INDEX_BASE_ZERO,
      shipdate_lt_rows, shipdate_lt_columns, shipdate_lt_IA, shipdate_lt_IA+1, shipdate_lt_JA, shipdate_lt_csr_values );

  ////////////////////////
  // START TIME MEASUREMENT
  ////////////////////////
  GET_TIME(global_time_start);

  // compute selection = shipdate_gt * shipdate_lt 
  sparse_status_t selection_result;

  selection_result = mkl_sparse_spmm ( SPARSE_OPERATION_NON_TRANSPOSE,
      shipdate_gt_matrix,
      shipdate_lt_matrix,
      &selection_matrix);

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

  // compute aggregation = quantity * bang
  sparse_status_t aggregation_result;
  struct matrix_descr descrA;
  descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

  bang_vector = (float*) mkl_malloc ((quantity_columns * sizeof(float)), MEM_LINE_SIZE );
  aggregation_vector = (float*) mkl_malloc ((quantity_columns * sizeof(float)), MEM_LINE_SIZE );

  aggregation_result = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, quantity_matrix , descrA, bang_vector, 1.0,  aggregation_vector);

  // compute intermediate_result = projection * selection
  sparse_status_t intermediate_result;

  intermediate_result = mkl_sparse_spmm ( SPARSE_OPERATION_NON_TRANSPOSE, 
      projection_matrix,
      selection_matrix, 
      &intermediate_matrix);

  // compute final_result = intermediate_result * aggregation
  sparse_status_t final_result;

  final_result = mkl_sparse_spmm ( SPARSE_OPERATION_NON_TRANSPOSE, 
      intermediate_matrix,
      quantity_matrix,
      &final_matrix);

  ////////////////////////
  // STOP TIME MEASUREMENT
  ////////////////////////
  GET_TIME(global_time_stop);
  writeResults( );

  sparse_index_base_t    indexing;
  mkl_sparse_s_export_csr ( final_matrix, &indexing, &final_rows, &final_columns, &final_IA, &final_IA_END, &final_JA, &final_csr_values);

  return 0;

}

