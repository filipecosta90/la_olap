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

  float* shipdate_gt_csr_values = NULL;
  MKL_INT* shipdate_gt_JA;
  MKL_INT* shipdate_gt_IA;
  MKL_INT shipdate_gt_rows;
  MKL_INT shipdate_gt_columns;
  MKL_INT shipdate_gt_nnz;

  float* shipdate_lt_csr_values = NULL;
  MKL_INT* shipdate_lt_JA;
  MKL_INT* shipdate_lt_IA;
  MKL_INT shipdate_lt_rows;
  MKL_INT shipdate_lt_columns;
  MKL_INT shipdate_lt_nnz;


  //read return flag
  tbl_read( "__tbl/lineitem.tbl" , 9, &returnFlag_nnz, &returnFlag_rows, &returnFlag_columns , &returnFlag_csr_values, &returnFlag_JA, &returnFlag_IA);

  //read line status
  tbl_read( "__tbl/lineitem.tbl" , 10, &lineStatus_nnz, &lineStatus_rows, &lineStatus_columns , &lineStatus_csr_values, &lineStatus_JA, &lineStatus_IA);

  //read quantity
  tbl_read( "__tbl/lineitem.tbl" , 5, &quantity_nnz, &quantity_rows, &quantity_columns , &quantity_csr_values, &quantity_JA, &quantity_IA);

  //read shipdate gt
  tbl_read_filter( "__tbl/lineitem.tbl" , 11, GREATER_EQ , "1998-08-28", &shipdate_gt_nnz, &shipdate_gt_rows, &shipdate_gt_columns , &shipdate_gt_csr_values, &shipdate_gt_JA, &shipdate_gt_IA);

  //read shipdate lt
  tbl_read_filter( "__tbl/lineitem.tbl" , 11, LESS_EQ , "1998-12-01", &shipdate_lt_nnz, &shipdate_lt_rows, &shipdate_lt_columns , &shipdate_lt_csr_values, &shipdate_lt_JA, &shipdate_lt_IA);

  float* selection_csr_values = NULL;
  MKL_INT* selection_JA;
  MKL_INT* selection_IA;
  MKL_INT selection_rows;
  MKL_INT selection_columns;
  MKL_INT selection_nnz;

  float* projection_csr_values = NULL;
  MKL_INT* projection_JA;
  MKL_INT* projection_IA;
  MKL_INT projection_rows;
  MKL_INT projection_columns;
  MKL_INT projection_nnz;

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

  float* final_csr_values = NULL;
  MKL_INT* final_JA;
  MKL_INT* final_IA;
  MKL_INT final_rows;
  MKL_INT final_columns;
  MKL_INT final_nnz;

  // compute C = A krao B
  /*  csr_kron( 
      A_csr_values, A_JA, A_IA, A_nnz, A_rows, A_columns, 
      B_csr_values, B_JA, B_IA, B_nnz, A_rows, B_columns, 
      &C_csr_values, &C_JA, &C_IA, &C_nnz, &C_rows, &C_columns 
      );
      */  
  return 0;
}
