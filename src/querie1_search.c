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
double global_time_selection,  global_time_selection_quantity, global_time_selection_aggregation , global_time_projection ;

void writeResults ( char* dataset ) {
	double total_time, selection_time, selection_quantity_time, selection_aggregation_time, projection_time, projection_selection_aggregation_time, final_time;
	selection_time = global_time_selection - global_time_start;
	selection_quantity_time = global_time_selection_quantity - global_time_selection;
	selection_aggregation_time = global_time_selection_aggregation - global_time_selection_quantity;
	projection_time = global_time_projection - global_time_selection_aggregation;
	projection_selection_aggregation_time = global_time_stop - global_time_projection;
	total_time = global_time_stop - global_time_start;
	char file_write[80];
	strcpy(file_write, "timing/timings_vec_");
	strcat(file_write, dataset);
	strcat(file_write, ".csv");

	FILE* stream = fopen(file_write, "a+");
	fprintf(stream, "%s, %lf, %lf, %lf, %lf, %lf, %lf\n",dataset, selection_time, selection_quantity_time, selection_aggregation_time, projection_time, projection_selection_aggregation_time , total_time);
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
	//////////////////////////////////////////
	//        CONVERT from CSR to CSC
	//////////////////////////////////////////

	MKL_INT job_csr_csc[8];
	// If job[0]=0, the matrix in the CSR format is converted to the CSC format;
	job_csr_csc[0] = 0;

	// job[1]
	// If job[1]=0, zero-based indexing for the matrix in CSR format is used;
	// if job[1]=1, one-based indexing for the matrix in CSR format is used.
	job_csr_csc[1] = 0;

	// job[2]
	// If job[2]=0, zero-based indexing for the matrix in the CSC format is used;
	// if job[2]=1, one-based indexing for the matrix in the CSC format is used.
	job_csr_csc[2] = 0;

	// job[5] - job indicator.
	// If job[5]=0, only arrays ja1, ia1 are filled in for the output storage.
	// If job[5]≠0, all output arrays acsc, ja1, and ia1 are filled in for the output storage.
	job_csr_csc[5] = 1;
	sparse_status_t status_convert_to_csc;

	/** ---------------------------------------------------------------------------
	 ** Return Flag Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	float* return_flag_csr_values = NULL;
	MKL_INT* return_flag_JA;
	MKL_INT* return_flag_IA;
	//CSC
	float* return_flag_csc_values = NULL;
	MKL_INT* return_flag_JA_csc;
	MKL_INT* return_flag_IA_csc;
	//COMMON
	MKL_INT return_flag_rows;
	MKL_INT return_flag_columns;
	MKL_INT return_flag_nnz;

	/** ---------------------------------------------------------------------------
	 ** Line Status Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	float* line_status_csr_values = NULL;
	MKL_INT* line_status_JA;
	MKL_INT* line_status_IA;
	//CSC
	float* line_status_csc_values = NULL;
	MKL_INT* line_status_JA_csc;
	MKL_INT* line_status_IA_csc;
	//COMMON
	MKL_INT line_status_rows;
	MKL_INT line_status_columns;
	MKL_INT line_status_nnz;

	/** ---------------------------------------------------------------------------
	 ** Quantity Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	float* quantity_csr_values = NULL;
	MKL_INT* quantity_JA;
	MKL_INT* quantity_IA;
	//CSC
	float* quantity_csc_values = NULL;
	MKL_INT* quantity_JA_csc;
	MKL_INT* quantity_IA_csc;
	//COMMON
	MKL_INT quantity_rows;
	MKL_INT quantity_columns;
	MKL_INT quantity_nnz;
	sparse_matrix_t  quantity_matrix;

	/** ---------------------------------------------------------------------------
	 ** Shipdate Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	float* shipdate_csr_values = NULL;
	MKL_INT* shipdate_JA;
	MKL_INT* shipdate_IA;
	//CSC
	float* shipdate_csc_values = NULL;
	MKL_INT* shipdate_JA_csc;
	MKL_INT* shipdate_IA_csc;
	//COMMON
	MKL_INT shipdate_rows;
	MKL_INT shipdate_columns;
	MKL_INT shipdate_nnz;
	sparse_matrix_t  shipdate_matrix;

	/* ---------------------------------------------------------------------------
	 ** Projection Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	float* projection_csr_values = NULL;
	MKL_INT* projection_JA;
	MKL_INT* projection_IA;
	//CSC
	float* projection_csc_values = NULL;
	MKL_INT* projection_JA_csc;
	MKL_INT* projection_IA_csc;

	//COMMON
	MKL_INT projection_rows;
	MKL_INT projection_columns;
	MKL_INT projection_nnz;
	sparse_matrix_t  projection_matrix;

	/* ---------------------------------------------------------------------------
	 ** Selection Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	float* selection_csr_values = NULL;
	MKL_INT* selection_JA;
	MKL_INT* selection_IA;
	//COMMON
	MKL_INT selection_rows;
	MKL_INT selection_columns;
	MKL_INT selection_nnz;
	sparse_matrix_t  selection_matrix;

	/* ---------------------------------------------------------------------------
	 ** Aggregation Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	float* aggregation_csr_values = NULL;
	MKL_INT* aggregation_JA;
	MKL_INT* aggregation_IA;
	//COMMON
	MKL_INT aggregation_rows;

	MKL_INT aggregation_nnz;

	/* ---------------------------------------------------------------------------
	 ** Intermediate Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	float* intermediate_csr_values = NULL;
	MKL_INT* intermediate_JA;
	MKL_INT* intermediate_IA;
	//COMMON
	MKL_INT intermediate_rows;
	MKL_INT intermediate_columns;
	MKL_INT intermediate_nnz;
	sparse_matrix_t  intermediate_matrix;

	/* ---------------------------------------------------------------------------
	 ** Vectors
	 ** -------------------------------------------------------------------------*/
	float* bang_vector;
	float* aggregation_vector;
	MKL_INT aggregation_vector_rows;


	float* intermediate_vector;
	MKL_INT intermediate_vector_rows;
	float* final_vector;


	//conversion status from csr arrays into mkl sparse_matrix_t 
	sparse_status_t status_to_csr;
	MKL_INT conversion_info = 0;
	/** ---------------------------------------------------------------------------
	 ** Populate Return Flag Matrix
	 ** -------------------------------------------------------------------------*/
	//read return flag
	tbl_read(
			table_file , 9, 
			&return_flag_nnz, &return_flag_rows, &return_flag_columns, 
			&return_flag_csr_values, &return_flag_JA, &return_flag_IA
		);

#ifdef D_DEBUGGING
	printf("return flag %d %d -- nnz: %d\n", return_flag_rows, return_flag_columns, return_flag_nnz );
#endif

	// Memory Allocation
	return_flag_csc_values = (float*) malloc ( return_flag_nnz * sizeof(float) );
	return_flag_JA_csc = (MKL_INT*) malloc ( return_flag_nnz * sizeof(MKL_INT) );
	return_flag_IA_csc = (MKL_INT*) malloc (( return_flag_nnz+1) * sizeof(MKL_INT));

	// Convert from CSR to CSC
	mkl_scsrcsc(job_csr_csc, &return_flag_nnz, return_flag_csr_values, return_flag_JA, return_flag_IA, return_flag_csc_values, return_flag_JA_csc, return_flag_IA_csc, &conversion_info);

#ifdef D_DEBUGGING
	csc_tbl_write(
			"return_flag_csc.txt",
			return_flag_csc_values, return_flag_JA_csc, return_flag_IA_csc,
			return_flag_nnz, return_flag_rows, return_flag_columns
		     );
#endif

	/** ---------------------------------------------------------------------------
	 ** Populate Line Status Matrix
	 ** -------------------------------------------------------------------------*/
	//read line status

	tbl_read(  
			table_file , 10, 
			&line_status_nnz, &line_status_rows, &line_status_columns , 
			&line_status_csr_values, &line_status_JA, &line_status_IA
		);

	// Memory Allocation
	line_status_csc_values = (float*) malloc ( line_status_nnz * sizeof(float) );
	line_status_JA_csc = (MKL_INT*) malloc ( line_status_nnz * sizeof(MKL_INT) );
	line_status_IA_csc = (MKL_INT*) malloc (( line_status_nnz+1) * sizeof(MKL_INT));

	// Convert from CSR to CSC
	mkl_scsrcsc(job_csr_csc, &line_status_nnz, line_status_csr_values, line_status_JA, line_status_IA, line_status_csc_values, line_status_JA_csc, line_status_IA_csc, &conversion_info);
#ifdef D_DEBUGGING
	printf("line status %d %d -- nnz: %d\n", line_status_rows, line_status_columns, line_status_nnz );

	csc_tbl_write(
			"line_status_csc.txt",
			line_status_csc_values, line_status_JA_csc, line_status_IA_csc,
			line_status_nnz, line_status_rows, line_status_columns
		     );
#endif

	/** ---------------------------------------------------------------------------
	 ** Populate Quantity Matrix
	 ** -------------------------------------------------------------------------*/
	//read quantity

	// measure
	tbl_read_measure(
			table_file , 5,
			&quantity_nnz,  &quantity_rows, &quantity_columns , 
			&quantity_csr_values, &quantity_JA, &quantity_IA
			);
#ifdef D_DEBUGGING
	printf("going to reshape quantity table\n");
#endif
	/*	csr_csr_square_reshape (
		&quantity_csr_values, &quantity_JA, &quantity_IA,
		&quantity_nnz, &quantity_rows, &quantity_columns,
		line_status_columns    
		);
		*/
#ifdef D_DEBUGGING
	csr_tbl_write(
			"quantity_test.txt",
			quantity_csr_values, quantity_JA, quantity_IA,
			quantity_nnz, quantity_rows, quantity_columns
		     );
#endif

	// Memory Allocation
	quantity_csc_values = (float*) malloc ( quantity_nnz * sizeof(float) );
	quantity_JA_csc = (MKL_INT*) malloc ( quantity_nnz * sizeof(MKL_INT) );
	quantity_IA_csc = (MKL_INT*) malloc (( quantity_nnz+1) * sizeof(MKL_INT));

	intermediate_csr_values = (float*) malloc ( (quantity_nnz) * sizeof(float) );
	intermediate_JA = (MKL_INT*) malloc ((quantity_nnz) * sizeof(MKL_INT) );
	intermediate_IA = (MKL_INT*) malloc (( quantity_nnz+1) * sizeof(MKL_INT));

	/** ---------------------------------------------------------------------------
	 ** Populate Shipdate Matrix
	 ** -------------------------------------------------------------------------*/
	//read shipdate
	tbl_read(
			table_file , 11,
			&shipdate_nnz, &shipdate_rows, &shipdate_columns ,
			&shipdate_csr_values, &shipdate_JA, &shipdate_IA
		);

#ifdef D_DEBUGGING
	printf("going to reshape quantity table\n");
#endif
	/*	csr_csr_square_reshape (
		&shipdate_csr_values, &shipdate_JA, &shipdate_IA,
		&shipdate_nnz, &shipdate_rows, &shipdate_columns,
		line_status_columns    
		);
		*/
#ifdef D_DEBUGGING
	csr_tbl_write(
			"shipdate_test_csr.txt",
			shipdate_csr_values, shipdate_JA, shipdate_IA,
			shipdate_nnz, shipdate_rows, shipdate_columns
		     );
#endif

	// Memory Allocation
	shipdate_csc_values = (float*) malloc ( shipdate_nnz * sizeof(float));
	shipdate_JA_csc = (MKL_INT*) malloc ( shipdate_nnz * sizeof(MKL_INT) );
	shipdate_IA_csc = (MKL_INT*) malloc ((shipdate_nnz+1) * sizeof(MKL_INT));

	// Convert from CSR to CSC
	mkl_scsrcsc(job_csr_csc, &shipdate_nnz, shipdate_csr_values, shipdate_JA, shipdate_IA, shipdate_csc_values, shipdate_JA_csc, shipdate_IA_csc, &conversion_info);

#ifdef D_DEBUGGING
	csc_tbl_write(
			"shipdate_test_csc.txt",
			shipdate_csc_values, shipdate_JA_csc, shipdate_IA_csc,
			shipdate_nnz, shipdate_rows, shipdate_columns
		     );
#endif

	//        convert via sparseBLAS API to Handle containing internal data for
	//        subsequent Inspector-executor Sparse BLAS operations.
	status_to_csr = mkl_sparse_s_create_csr ( 
			&shipdate_matrix , SPARSE_INDEX_BASE_ZERO,      shipdate_rows, shipdate_columns, 
			shipdate_IA, shipdate_IA+1, shipdate_JA, shipdate_csr_values
			);

	// Convert from CSR to CSC
	mkl_scsrcsc(job_csr_csc, &quantity_nnz, quantity_csr_values, quantity_JA, quantity_IA, quantity_csc_values, quantity_JA_csc, quantity_IA_csc, &conversion_info);

#ifdef D_DEBUGGING
	csc_tbl_write(
			"quantity_test_csc.txt",
			quantity_csc_values, quantity_JA_csc, quantity_IA_csc,
			quantity_nnz, quantity_rows, quantity_columns
		     );
#endif

	aggregation_vector_rows = quantity_columns;
	aggregation_vector = (float*) malloc ((aggregation_vector_rows+1) * sizeof(float));

	//        convert via sparseBLAS API to Handle containing internal data for 
	//        subsequent Inspector-executor Sparse BLAS operations.
	status_to_csr = mkl_sparse_s_create_csr ( &quantity_matrix , SPARSE_INDEX_BASE_ZERO, 
			quantity_rows, quantity_columns, quantity_IA, quantity_IA+1, quantity_JA, quantity_csr_values );
#ifdef D_DEBUGGING
	check_errors(status_to_csr);
#endif

	/** ---------------------------------------------------------------------------
	 ** Auxiliar Vars
	 ** -------------------------------------------------------------------------*/

	sparse_status_t selection_result;
	sparse_status_t aggregation_result;
	sparse_status_t intermediate_result;
	sparse_status_t final_result;
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

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
	 ** ---------------------------------------------------------------------------
	 ** -------------------------------------------------------------------------*/

		mkl_sparse_optimize(quantity_matrix);
/*
		MKL_INT max_threads;
		max_threads = mkl_get_max_threads();
		mkl_set_num_threads(max_threads);
		*/
	//	printf("** Setted max thread on MKL to %d\n", max_threads);
#ifdef D_DEBUGGING
	printf("** START TIME MEASUREMENT\n");
#endif
//	sleep(5);
	GET_TIME(global_time_start);

	csc_to_csr_mx_selection_and(
			shipdate_csc_values, shipdate_JA_csc, shipdate_IA_csc,
			shipdate_nnz, shipdate_rows, shipdate_columns,
			GREATER_EQ , "1998-08-28", LESS_EQ , "1998-12-01",
			&selection_csr_values, &selection_JA, &selection_IA,
			&selection_nnz, &selection_rows, &selection_columns
			);

	intermediate_vector_rows = selection_columns;
printf("#######selecao\n");
	print_csr(
			selection_csr_values, selection_JA, selection_IA,
			selection_nnz, selection_rows, selection_columns
		 );


	intermediate_vector = (float*) malloc ( (intermediate_vector_rows+1) * sizeof(float));
	bang_vector = (float*) malloc ( (quantity_columns+1) * sizeof(float));
	for (int pos =0; pos < line_status_columns ; pos++){
		bang_vector[pos] = 1.0; 
	}
#ifdef D_DEBUGGING
	csr_tbl_write(
			"selection_test.txt",
			selection_csr_values, selection_JA, selection_IA,
			selection_nnz, selection_rows, selection_columns
		     );
#endif

	status_to_csr = mkl_sparse_s_create_csr ( &selection_matrix , SPARSE_INDEX_BASE_ZERO, selection_rows, selection_columns, selection_IA, selection_IA+1, selection_JA, selection_csr_values );
#ifdef D_DEBUGGING
	printf("quantity rows %d columns %d nnz %d\n", quantity_rows, quantity_columns, quantity_nnz);
	printf("selection rows %d columns %d nnz %d\n", selection_rows, selection_columns, selection_nnz);
	check_errors(status_to_csr);
#endif
	GET_TIME(global_time_selection);

	// compute intermediate = selection * aggregation
	intermediate_result = mkl_sparse_spmm ( SPARSE_OPERATION_NON_TRANSPOSE,
			selection_matrix,
			quantity_matrix,
			&intermediate_matrix);

#ifdef D_DEBUGGING
	check_errors(intermediate_result);
	printf(" // compute aggregation = quantity * bang\n");
#endif
	// mkl_sparse_optimize(projection_matrix);
	GET_TIME(global_time_selection_quantity);

	// compute aggregation = quantity * bang
	// results in a vector


	aggregation_result = mkl_sparse_s_mv (
			SPARSE_OPERATION_NON_TRANSPOSE, 1.0, intermediate_matrix , descrA, bang_vector, 0.0,  intermediate_vector
			);

	GET_TIME(global_time_selection_aggregation);

#ifdef D_DEBUGGING
	csr_measure_vector_write(
			"intermediate_test.txt",
			intermediate_vector, intermediate_vector_rows
			);
#endif
	// compute projection = return_flag krao line_status

#ifdef D_DEBUGGING

	print_csr(
			return_flag_csr_values, return_flag_JA, return_flag_IA,
			return_flag_nnz, return_flag_rows, return_flag_columns
		 );

	print_csr(
			line_status_csr_values, line_status_JA, line_status_IA,
			line_status_nnz, line_status_rows, line_status_columns
		 );

	print_csc(
			projection_csc_values, projection_JA_csc, projection_IA_csc,
			projection_nnz, projection_rows, projection_columns
		 );
printf("INTERESSA!!!!\n");
	print_csr(
			projection_csr_values, projection_JA, projection_IA,
			projection_nnz, projection_rows, projection_columns
		 );


	csc_tbl_write(
			"projection_test_csc.txt",
			projection_csc_values, projection_JA_csc, projection_IA_csc,
			projection_nnz, projection_rows, projection_columns
		     );

	csr_tbl_write(
			"projection_test_csr.txt",
			projection_csr_values, projection_JA, projection_IA,
			projection_nnz, projection_rows, projection_columns
		     );
	printf("projection nnz: %d\n", projection_nnz);
	printf("projection rows: %d\n", projection_rows);
	printf("projection columns: %d\n", projection_columns);
#endif
	status_to_csr = mkl_sparse_s_create_csr ( &projection_matrix , SPARSE_INDEX_BASE_ZERO, projection_rows, projection_columns, projection_IA, projection_IA+1, projection_JA, projection_csr_values );

	GET_TIME(global_time_projection);

	// compute final_result = intermediate_result * aggregation

	final_result = mkl_sparse_s_mv (
			SPARSE_OPERATION_NON_TRANSPOSE, 1.0, projection_matrix , descrA, intermediate_vector, 0.0,  final_vector
			);
#ifdef D_DEBUGGING
	printf("end of compute\n");
#endif
	////////////////////////
	// STOP TIME MEASUREMENT
	////////////////////////
	GET_TIME(global_time_stop);
#ifdef D_DEBUGGING
	printf("** STOP TIME MEASUREMENT\n");
#endif
	writeResults( argv[1] );

	csr_measure_vector_write(
			"final_test.txt",
			final_vector, projection_rows
			);

	return 0;
}

