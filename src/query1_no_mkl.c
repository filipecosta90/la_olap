/* ---------------------------------------------------------------------------
 **    Filename: querie1_search.c
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
 ** University of Minho, High Performance Computing Dpt. , May 2016
 ** -------------------------------------------------------------------------*/

#include <stdio.h>
#include <glib.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mkl_spblas.h"
#include "mkl_types.h"
#include "mkl.h"
#include "olap_search.h"
#include "timer.h"
//sleep
#include <unistd.h>


double global_time_start, global_time_stop;
// intermediate timers
double global_time_projection, global_time_selection,  global_time_projection_selection, global_time_projection_selection_aggregation , global_time_bang;

void writeResults ( char* dataset ) {
  double total_time, projection_time, selection_time, projection_selection_time, projection_selection_aggregation_time, bang_time, final_time;

  //1st opp
  projection_time = global_time_projection - global_time_start;

  //2nd opp
  selection_time = global_time_selection - global_time_projection;

  //3rd opp
  projection_selection_time = global_time_projection_selection - global_time_selection;

  //4th opp
  projection_selection_aggregation_time = global_time_projection_selection_aggregation - global_time_projection_selection;

  //5th opp
  bang_time = global_time_bang - global_time_projection_selection_aggregation;

  total_time = global_time_stop - global_time_start;
  char file_write[80];
  strcpy(file_write, "timing/timings_vec_nomkl_");
  strcat(file_write, dataset);
  strcat(file_write, ".csv");

  FILE* stream = fopen(file_write, "a+");
  fprintf(stream, "%s, %lf, %lf, %lf, %lf, %lf, %lf\n",dataset, projection_time, selection_time, projection_selection_time, projection_selection_aggregation_time, bang_time , total_time);
  fclose(stream);
}

int main( int argc, char* argv[]){

  char table_file[80];
  strcpy(table_file, "__tbl/lineitem_");
  strcat(table_file, argv[1]);
  strcat(table_file, ".tbl");

#ifdef D_DEBUGGING
  printf("going to read results from %s\n", table_file);
#endif


  /** ---------------------------------------------------------------------------
   ** Return Flag Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* return_flag_csc_values = NULL;
  int* return_flag_row_ind;
  int* return_flag_col_ptr;
  //COMMON
  int return_flag_n_rows;
  int return_flag_n_cols;
  int return_flag_n_nnz;

  /** ---------------------------------------------------------------------------
   ** Line Status Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* line_status_csc_values = NULL;
  int* line_status_row_ind;
  int* line_status_col_ptr;
  //COMMON
  int line_status_n_rows;
  int line_status_n_cols;
  int line_status_n_nnz;

  /** ---------------------------------------------------------------------------
   ** Quantity Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* quantity_csc_values = NULL;
  int* quantity_row_ind;
  int* quantity_col_ptr;
  //COMMON
  int quantity_n_rows;
  int quantity_n_cols;
  int quantity_n_nnz;

  /** ---------------------------------------------------------------------------
   ** Shipdate Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* shipdate_csc_values = NULL;
  int* shipdate_row_ind;
  int* shipdate_col_ptr;
  //COMMON
  int shipdate_n_rows;
  int shipdate_n_cols;
  int shipdate_n_nnz;

  /* ---------------------------------------------------------------------------
   ** Projection Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* projection_csc_values = NULL;
  int* projection_row_ind;
  int* projection_col_ptr;

  //COMMON
  int projection_n_rows;
  int projection_n_cols;
  int projection_n_nnz;

  /* ---------------------------------------------------------------------------
   ** Selection Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* selection_csc_values = NULL;
  int* selection_row_ind;
  int* selection_col_ptr;
  //COMMON
  int selection_n_rows;
  int selection_n_cols;
  int selection_n_nnz;

  /* ---------------------------------------------------------------------------
   ** Projection Selection Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* projection_selection_csc_values;
  int* projection_selection_row_ind;
  int* projection_selection_col_ind;
  //COMMON
  int  projection_selection_n_nnz;
  int projection_selection_n_rows;
  int  projection_selection_n_columns;

  /* ---------------------------------------------------------------------------
   ** Intermediate Matrix
   ** -------------------------------------------------------------------------*/
  //CSR
  float* intermediate_csc_values = NULL;
  int* intermediate_JA;
  int* intermediate_IA;
  //COMMON
  int intermediate_n_rows;
  int intermediate_n_cols;
  int intermediate_n_nnz;
  sparse_matrix_t  intermediate_matrix;

  /* ---------------------------------------------------------------------------
   ** Vectors
   ** -------------------------------------------------------------------------*/
  float* bang_vector;
  float* aggregation_vector;
  int aggregation_vector_n_rows;


  float* intermediate_vector;
  int intermediate_vector_n_rows;

  float* final_vector;

  int number_elements = tbl_get_number_elements (table_file);

#ifdef D_DEBUGGING
  printf("tbl file has %d elements\n", number_elements);
#endif
  /** ---------------------------------------------------------------------------
   ** =========================== END OF DECLARATIONS ===========================
   ** -------------------------------------------------------------------------*/

  /** ---------------------------------------------------------------------------
   ** Populate Return Flag Matrix
   ** -------------------------------------------------------------------------*/
  //read return flag
  //bitmap matrix
  tbl_read_csc(
      table_file , 9, number_elements,
      &return_flag_n_nnz, &return_flag_n_rows, &return_flag_n_cols, 
      &return_flag_csc_values, &return_flag_row_ind, &return_flag_col_ptr
      );

  /** ---------------------------------------------------------------------------
   ** Populate Line Status Matrix
   ** -------------------------------------------------------------------------*/
  //read line status
  //bitmap matrix
  tbl_read_csc(
      table_file , 10, number_elements,
      &line_status_n_nnz, &line_status_n_rows, &line_status_n_cols , 
      &line_status_csc_values, &line_status_row_ind, &line_status_col_ptr
      );

  /** ---------------------------------------------------------------------------
   ** Populate Quantity Matrix
   ** -------------------------------------------------------------------------*/
  //read quantity
  // measure
  tbl_read_csc_measure(
      table_file , 5, number_elements,
      &quantity_n_nnz,  &quantity_n_rows, &quantity_n_cols , 
      &quantity_csc_values, &quantity_row_ind, &quantity_col_ptr
      );

  /** ---------------------------------------------------------------------------
   ** Populate Shipdate Matrix
   ** -------------------------------------------------------------------------*/
  //read shipdate
  //bitmap matrix
  tbl_read_csc(
      table_file , 11, number_elements,
      &shipdate_n_nnz, &shipdate_n_rows, &shipdate_n_cols ,
      &shipdate_csc_values, &shipdate_row_ind, &shipdate_col_ptr
      );


  aggregation_vector_n_rows = quantity_n_cols;
  aggregation_vector = (float*) _mm_malloc ((aggregation_vector_n_rows+1) * sizeof(float), MEM_LINE_SIZE );


  /** ---------------------------------------------------------------------------
   ** Populate Vectors
   ** -------------------------------------------------------------------------*/

  final_vector = (float*) malloc ( (quantity_n_cols+1) * sizeof(float));

  /** ---------------------------------------------------------------------------
   ** ---------------------------------------------------------------------------
   ** ---------------------------------------------------------------------------
   ** ---------------------------------------------------------------------------
   ** ---------------------------------------------------------------------------
   ** START TIME MEASUREMENT
   ** ---------------------------------------------------------------------------
   ** ---------------------------------------------------------------------------
   ** ---------------------------------------------------------------------------
   ** ---------------------------------------------------------------------------
   ** -------------------------------------------------------------------------*/


#ifdef D_DEBUGGING
  printf("** START TIME MEASUREMENT\n");
  //sleep(5);
#endif

  GET_TIME(global_time_start);


  csc_csc_krao(
      return_flag_csc_values, return_flag_row_ind, return_flag_col_ptr,
      return_flag_n_nnz, return_flag_n_rows, return_flag_n_cols,

      line_status_csc_values, line_status_row_ind, line_status_col_ptr ,
      line_status_n_nnz, line_status_n_rows, line_status_n_cols,

      &projection_csc_values, &projection_row_ind, &projection_col_ptr,
      &projection_n_nnz, &projection_n_rows, &projection_n_cols
      );


  GET_TIME(global_time_projection);

  csc_to_csc_mx_selection_and(
      shipdate_csc_values, shipdate_row_ind, shipdate_col_ptr,
      shipdate_n_nnz, shipdate_n_rows, shipdate_n_cols,
      GREATER_EQ , "1998-08-28", LESS_EQ , "1998-12-01",
      &selection_csc_values, &selection_row_ind, &selection_col_ptr,
      &selection_n_nnz, &selection_n_rows, &selection_n_cols
      );

#ifdef D_DEBUGGING
  printf("selection rows %d columns %d nnz %d\n", selection_n_rows, selection_n_cols, selection_n_nnz);

  print_csc (
      selection_csc_values, selection_row_ind, selection_col_ptr,
      selection_n_nnz, selection_n_rows, selection_n_cols
      );
#endif

  GET_TIME(global_time_selection);

  csc_csc_mm(
      projection_csc_values, projection_row_ind, projection_col_ind,
      projection_n_nnz, projection_n_rows, projection_n_columns,

      selection_csc_values, selection_row_ind, selection_col_ptr,
      selection_n_nnz, selection_n_rows, selection_n_cols,

      &projection_selection_csc_values, &projection_selection_row_ind, &projection_selection_col_ind,
      &projection_selection_n_nnz, &projection_selection_n_rows, &projection_selection_n_columns
      );

  GET_TIME(global_time_projection_selection);

  GET_TIME(global_time_projection_selection_aggregation);

  GET_TIME(global_time_bang);

  ////////////////////////
  // STOP TIME MEASUREMENT
  ////////////////////////
  GET_TIME(global_time_stop);

#ifdef D_DEBUGGING
  printf("** STOP TIME MEASUREMENT\n");
#endif
  writeResults( argv[1] );

  return 0;
}

