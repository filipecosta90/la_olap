/* ---------------------------------------------------------------------------
 **    Filename: querie1_data.c
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
  /*

     printf("reading return flag\n");
     tbl_read( "__tbl_16/lineitem.tbl" , 9, &returnFlag_nnz, &returnFlag_rows, &returnFlag_columns, &returnFlag_csr_values, &returnFlag_JA, &returnFlag_IA);
     printf("generating return flag\n");
     convert_and_write_to_mx("__quark_mx_16/return_flag_16.mx", returnFlag_csr_values, returnFlag_JA, returnFlag_IA, returnFlag_nnz, returnFlag_rows, returnFlag_columns);
     mkl_free(returnFlag_csr_values);
     mkl_free(returnFlag_JA);
     mkl_free(returnFlag_IA);

  //read quantity
  printf("reading quantity\n");
  tbl_read_measure( "__tbl_16/lineitem.tbl" , 5, &quantity_nnz, &quantity_rows, &quantity_columns , &quantity_csr_values, &quantity_JA, &quantity_IA);
  printf("generating quantity\n");
  convert_and_write_to_mx("__quark_mx_16/quantity_16.mx", quantity_csr_values, quantity_JA, quantity_IA, quantity_nnz, quantity_rows, quantity_columns);
  mkl_free(quantity_csr_values);
  mkl_free(quantity_JA);
  mkl_free(quantity_IA);



  //read line status

  printf("reading line status\n");
  tbl_read( "__tbl_16/lineitem.tbl" , 10, &lineStatus_nnz, &lineStatus_rows, &lineStatus_columns , &lineStatus_csr_values, &lineStatus_JA, &lineStatus_IA);
  printf("generating line status\n");
  convert_and_write_to_mx("__quark_mx_16/line_status_16.mx", lineStatus_csr_values, lineStatus_JA, lineStatus_IA, lineStatus_nnz, lineStatus_rows, lineStatus_columns);
  mkl_free(lineStatus_csr_values);
  mkl_free(lineStatus_JA);
  mkl_free(lineStatus_IA);



  //read shipdate gt
  printf("reading shipdate gt\n");
  tbl_read_filter( "__tbl_16/lineitem.tbl" , 11, GREATER_EQ , "1998-08-28",
  &shipdate_gt_nnz, &shipdate_gt_rows, &shipdate_gt_columns , &shipdate_gt_csr_values, &shipdate_gt_JA, &shipdate_gt_IA);
  printf("generating shipdate gt\n");
  convert_and_write_to_mx("__quark_mx_16/shipdate_gt_16.mx", shipdate_gt_csr_values, shipdate_gt_JA, shipdate_gt_IA, shipdate_gt_nnz, shipdate_gt_rows, shipdate_gt_columns);
  mkl_free(shipdate_gt_csr_values);
  mkl_free(shipdate_gt_JA);
  mkl_free(shipdate_gt_IA);



*/

  //read shipdate lt
  printf("reading shipdate lt\n");
  tbl_read_filter( "__tbl_16/lineitem.tbl" , 11, LESS_EQ , "1998-12-01",
      &shipdate_lt_nnz, &shipdate_lt_rows, &shipdate_lt_columns , &shipdate_lt_csr_values, &shipdate_lt_JA, &shipdate_lt_IA);
  printf("generating shipdate lt\n");
  convert_and_write_to_mx("__quark_mx_16/shipdate_lt_16.mx", shipdate_lt_csr_values, shipdate_lt_JA, shipdate_lt_IA, shipdate_lt_nnz, shipdate_lt_rows, shipdate_lt_columns);
  mkl_free(shipdate_lt_csr_values);
  mkl_free(shipdate_lt_JA);
  mkl_free(shipdate_lt_IA);
  return 0;

}

