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
  int* return_flag_JA_csc;
  int* return_flag_IA_csc;
  //COMMON
  int return_flag_rows;
  int return_flag_columns;
  int return_flag_nnz;

  /** ---------------------------------------------------------------------------
   ** Line Status Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* line_status_csc_values = NULL;
  int* line_status_JA_csc;
  int* line_status_IA_csc;
  //COMMON
  int line_status_rows;
  int line_status_columns;
  int line_status_nnz;

  /** ---------------------------------------------------------------------------
   ** Quantity Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* quantity_csc_values = NULL;
  int* quantity_JA_csc;
  int* quantity_IA_csc;
  //COMMON
  int quantity_rows;
  int quantity_columns;
  int quantity_nnz;

  /** ---------------------------------------------------------------------------
   ** Shipdate Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* shipdate_csc_values = NULL;
  int* shipdate_JA_csc;
  int* shipdate_IA_csc;
  //COMMON
  int shipdate_rows;
  int shipdate_columns;
  int shipdate_nnz;

  /* ---------------------------------------------------------------------------
   ** Projection Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* projection_csc_values = NULL;
  int* projection_JA_csc;
  int* projection_IA_csc;

  //COMMON
  int projection_rows;
  int projection_columns;
  int projection_nnz;

  /* ---------------------------------------------------------------------------
   ** Selection Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* selection_csc_values = NULL;
  int* selection_JA_csc;
  int* selection_IA_csc;
  //COMMON
  int selection_rows;
  int selection_columns;
  int selection_nnz;

  /* ---------------------------------------------------------------------------
   ** Aggregation Matrix
   ** -------------------------------------------------------------------------*/
  //CSC
  float* aggregation_csc_values = NULL;
  int* aggregation_JA_csc;
  int* aggregation_IA_csc;
  //COMMON
  int aggregation_rows;

  int aggregation_nnz;

  /* ---------------------------------------------------------------------------
   ** Intermediate Matrix
   ** -------------------------------------------------------------------------*/
  //CSR
  float* intermediate_csc_values = NULL;
  int* intermediate_JA;
  int* intermediate_IA;
  //COMMON
  int intermediate_rows;
  int intermediate_columns;
  int intermediate_nnz;
  sparse_matrix_t  intermediate_matrix;

  /* ---------------------------------------------------------------------------
   ** Vectors
   ** -------------------------------------------------------------------------*/
  float* bang_vector;
  float* aggregation_vector;
  int aggregation_vector_rows;


  float* intermediate_vector;
  int intermediate_vector_rows;

  float* final_vector;

  int number_elements = tbl_get_number_elements (table_file);

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
      &return_flag_nnz, &return_flag_rows, &return_flag_columns, 
      &return_flag_csc_values, &return_flag_JA_csc, &return_flag_IA_csc
      );

  /** ---------------------------------------------------------------------------
   ** Populate Line Status Matrix
   ** -------------------------------------------------------------------------*/
  //read line status
  //bitmap matrix
  tbl_read_csc(
      table_file , 10, number_elements,
      &line_status_nnz, &line_status_rows, &line_status_columns , 
      &line_status_csc_values, &line_status_JA_csc, &line_status_IA_csc
      );

  /** ---------------------------------------------------------------------------
   ** Populate Quantity Matrix
   ** -------------------------------------------------------------------------*/
  //read quantity
  // measure
  tbl_read_csc_measure(
      table_file , 5, number_elements,
      &quantity_nnz,  &quantity_rows, &quantity_columns , 
      &quantity_csc_values, &quantity_JA_csc, &quantity_IA_csc
      );

  /** ---------------------------------------------------------------------------
   ** Populate Shipdate Matrix
   ** -------------------------------------------------------------------------*/
  //read shipdate
  //bitmap matrix
  tbl_read(
      table_file , 11, number_elements,
      &shipdate_nnz, &shipdate_rows, &shipdate_columns ,
      &shipdate_csc_values, &shipdate_JA_csc, &shipdate_IA_csc
      );


  aggregation_vector_rows = quantity_columns;
  aggregation_vector = (float*) _mm_malloc ((aggregation_vector_rows+1) * sizeof(float), MEM_LINE_SIZE );


  /** ---------------------------------------------------------------------------
   ** Populate Vectors
   ** -------------------------------------------------------------------------*/

  final_vector = (float*) malloc ( (quantity_columns+1) * sizeof(float));

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



  GET_TIME(global_time_projection);


  csc_to_csc_mx_selection_and(
      shipdate_csc_values, shipdate_JA_csc, shipdate_IA_csc,
      shipdate_nnz, shipdate_rows, shipdate_columns,
      GREATER_EQ , "1998-08-28", LESS_EQ , "1998-12-01",
      &selection_csc_values, &selection_JA_csc, &selection_IA_csc,
      &selection_nnz, &selection_rows, &selection_columns
      );

#ifdef D_DEBUGGING
  printf("selection rows %d columns %d nnz %d\n", selection_rows, selection_columns, selection_nnz);

  print_csc (
      selection_csc_values, selection_JA_csc, selection_IA_csc,
      selection_nnz, selection_rows, selection_columns
      );
#endif

  GET_TIME(global_time_selection);

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

