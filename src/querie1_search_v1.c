/* ---------------------------------------------------------------------------
 **    Filename: querie1_search_v1.c
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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mkl_spblas.h"
#include "mkl_types.h"
#include "mkl.h"
#include "olap_search.h"
#include "timer.h"


double global_time_start, global_time_stop;
double global_time_selection, global_time_projection, global_time_aggregation, global_time_intermediate,global_time_projection_final, global_time_projection_lw, global_time_projection_rq;

void writeResults ( char* dataset ) {
	double total_time, selection_time, projection_time, aggregation_time, intermediate_time, final_time;

	selection_time = global_time_selection - global_time_start;
	projection_time = global_time_projection - global_time_start;
	aggregation_time = global_time_aggregation - global_time_start;
	intermediate_time = global_time_intermediate - global_time_start;
	total_time = global_time_stop - global_time_start;
	char file_write[80];
	strcpy(file_write, "timing/timings_vec_");
	strcat(file_write, dataset);
	strcat(file_write, ".csv");

	FILE* stream = fopen(file_write, "a+");
	fprintf(stream, "%s,%f,%f,%f,%f,%f\n",dataset, selection_time, projection_time, aggregation_time, intermediate_time, total_time);
	fclose(stream);
}



int main( int argc, char* argv[]){


	char table_file[80];
	strcpy(table_file, "__tbl/lineitem_");
	strcat(table_file, argv[1]);
	strcat(table_file, ".tbl");

	printf("going to read results from %s\n", table_file);

	MKL_INT* quark_start_end;
	quark_start_end = (MKL_INT*) mkl_malloc (( MEM_LINE_SIZE * sizeof(MKL_INT) ), MEM_LINE_SIZE );
	MKL_INT  quark_distinct_tables = 0;
	quark_start_end[0] = 0;


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
	__declspec(align(MEM_LINE_SIZE))  float* return_flag_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* return_flag_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* return_flag_IA;
	//CSC
	__declspec(align(MEM_LINE_SIZE))  float* return_flag_csc_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* return_flag_JA1_csc;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* return_flag_IA1_csc;
	//COMMON
	MKL_INT return_flag_rows;
	MKL_INT return_flag_columns;
	MKL_INT return_flag_nnz;

	/** ---------------------------------------------------------------------------
	 ** Line Status Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	__declspec(align(MEM_LINE_SIZE))  float* line_status_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* line_status_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* line_status_IA;
	//CSC
	__declspec(align(MEM_LINE_SIZE))  float* line_status_csc_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* line_status_JA1_csc;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* line_status_IA1_csc;
	//COMMON
	MKL_INT line_status_rows;
	MKL_INT line_status_columns;
	MKL_INT line_status_nnz;

	/** ---------------------------------------------------------------------------
	 ** Quantity Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	__declspec(align(MEM_LINE_SIZE))  float* quantity_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* quantity_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* quantity_IA;
	//CSC
	__declspec(align(MEM_LINE_SIZE))  float* quantity_csc_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* quantity_JA1_csc;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* quantity_IA1_csc;
	//COMMON
	MKL_INT quantity_rows;
	MKL_INT quantity_columns;
	MKL_INT quantity_nnz;
	sparse_matrix_t  quantity_matrix;

	/** ---------------------------------------------------------------------------
	 ** Shipdate Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	__declspec(align(MEM_LINE_SIZE))  float* shipdate_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* shipdate_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* shipdate_IA;
	//CSR
	__declspec(align(MEM_LINE_SIZE))  float* shipdate_csc_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* shipdate_JA1_csc;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* shipdate_IA1_csc;
	//COMMON
	MKL_INT shipdate_rows;
	MKL_INT shipdate_columns;
	MKL_INT shipdate_nnz;
	sparse_matrix_t  shipdate_matrix;

	/* ---------------------------------------------------------------------------
	 ** Projection return_flag quantity Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	__declspec(align(MEM_LINE_SIZE))  float* projection_rq_csc_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* projection_rq_JA1_csc;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* projection_rq_IA1_csc;
	//COMMON
	MKL_INT projection_rq_rows;
	MKL_INT projection_rq_columns;
	MKL_INT projection_rq_nnz;
	sparse_matrix_t  projection_rq_matrix;

	/* ---------------------------------------------------------------------------
	 ** Projection return_flag quantity Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	__declspec(align(MEM_LINE_SIZE))  float* projection_final_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* projection_final_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* projection_final_IA;
	//CSC
	__declspec(align(MEM_LINE_SIZE))  float* projection_final_csc_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* projection_final_JA1_csc;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* projection_final_IA1_csc;

	//COMMON
	MKL_INT projection_final_rows;
	MKL_INT projection_final_columns;
	MKL_INT projection_final_nnz;
	sparse_matrix_t  projection_final_matrix;

	/* ---------------------------------------------------------------------------
	 ** Projection line_status selection Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	__declspec(align(MEM_LINE_SIZE))  float* projection_lw_csc_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* projection_lw_JA1_csc;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* projection_lw_IA1_csc;
	//COMMON
	MKL_INT projection_lw_rows;
	MKL_INT projection_lw_columns;
	MKL_INT projection_lw_nnz;
	sparse_matrix_t  projection_lw_matrix;

	/* ---------------------------------------------------------------------------
	 ** Selection Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	__declspec(align(MEM_LINE_SIZE))  float* selection_csc_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* selection_JA1_csc;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* selection_IA1_csc;
	//COMMON
	MKL_INT selection_rows;
	MKL_INT selection_columns;
	MKL_INT selection_nnz;
	sparse_matrix_t  selection_matrix;

	/* ---------------------------------------------------------------------------
	 ** Aggregation Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	__declspec(align(MEM_LINE_SIZE)) float* aggregation_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* aggregation_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* aggregation_IA;
	//COMMON
	MKL_INT aggregation_rows;
	MKL_INT aggregation_columns;
	MKL_INT aggregation_nnz;

	/* ---------------------------------------------------------------------------
	 ** Intermediate Matrix
	 ** -------------------------------------------------------------------------*/
	//CSR
	__declspec(align(MEM_LINE_SIZE)) float* intermediate_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* intermediate_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* intermediate_IA;
	//COMMON
	MKL_INT intermediate_rows;
	MKL_INT intermediate_columns;
	MKL_INT intermediate_nnz;
	sparse_matrix_t  intermediate_matrix;

	/* ---------------------------------------------------------------------------
	 ** Vectors
	 ** -------------------------------------------------------------------------*/
	__declspec(align(MEM_LINE_SIZE))  float* bang_vector;
	__declspec(align(MEM_LINE_SIZE))  float* aggregation_vector;
	__declspec(align(MEM_LINE_SIZE))  float* intermediate_vector;
	__declspec(align(MEM_LINE_SIZE))  float* final_vector;


	//conversion status from csr arrays into mkl sparse_matrix_t
	sparse_status_t status_to_csr;

	/** ---------------------------------------------------------------------------
	 ** Populate Return Flag Matrix
	 ** -------------------------------------------------------------------------*/
	//read return flag
	tbl_read(
			table_file , 9,
			&return_flag_nnz, &return_flag_rows, &return_flag_columns,
			&return_flag_csr_values, &return_flag_JA, &return_flag_IA
		);

	// Memory Allocation
	return_flag_csc_values = (float*) mkl_malloc (( return_flag_nnz * sizeof(float) ), MEM_LINE_SIZE );
	return_flag_JA1_csc = (MKL_INT*) mkl_malloc (( return_flag_nnz * sizeof(MKL_INT) ), MEM_LINE_SIZE );
	return_flag_IA1_csc = (MKL_INT*) mkl_malloc ((( return_flag_nnz+1) * sizeof(MKL_INT)), MEM_LINE_SIZE );

	// Convert from CSR to CSC
	mkl_scsrcsc(job_csr_csc, &return_flag_nnz, return_flag_csr_values, return_flag_JA, return_flag_IA, return_flag_csc_values, return_flag_JA1_csc, return_flag_IA1_csc, &status_convert_to_csc);
	printf("conversion of return flag matrix from CSR to CSC ok?\n\t");
	check_errors(status_convert_to_csc);

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
	line_status_csc_values = (float*) mkl_malloc (( line_status_nnz * sizeof(float) ), MEM_LINE_SIZE );
	line_status_JA1_csc = (MKL_INT*) mkl_malloc (( line_status_nnz * sizeof(MKL_INT) ), MEM_LINE_SIZE );
	line_status_IA1_csc = (MKL_INT*) mkl_malloc ((( line_status_nnz+1) * sizeof(MKL_INT)), MEM_LINE_SIZE );

	// Convert from CSR to CSC
	mkl_scsrcsc(job_csr_csc, &line_status_nnz, line_status_csr_values, line_status_JA, line_status_IA, line_status_csc_values, line_status_JA1_csc, line_status_IA1_csc, &status_convert_to_csc);
	printf("conversion of linestatus matrix from CSR to CSC ok?\n\t");
	check_errors(status_convert_to_csc);

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

	// Memory Allocation
	quantity_csc_values = (float*) mkl_malloc (( quantity_nnz * sizeof(float) ), MEM_LINE_SIZE );
	quantity_JA1_csc = (MKL_INT*) mkl_malloc (( quantity_nnz * sizeof(MKL_INT) ), MEM_LINE_SIZE );
	quantity_IA1_csc = (MKL_INT*) mkl_malloc ((( quantity_nnz+1) * sizeof(MKL_INT)), MEM_LINE_SIZE );

	intermediate_csr_values = (float*) mkl_malloc (( (quantity_nnz) * sizeof(float) ), MEM_LINE_SIZE );
	intermediate_JA = (MKL_INT*) mkl_malloc (((quantity_nnz) * sizeof(MKL_INT) ), MEM_LINE_SIZE );
	intermediate_IA = (MKL_INT*) mkl_malloc ((( quantity_nnz+1) * sizeof(MKL_INT)), MEM_LINE_SIZE );

	// Convert from CSR to CSC
	mkl_scsrcsc(job_csr_csc, &quantity_nnz, quantity_csr_values, quantity_JA, quantity_IA, quantity_csc_values, quantity_JA1_csc, quantity_IA1_csc, &status_convert_to_csc);
	printf("conversion of quantity matrix from CSR to CSC ok?\n\t");
	check_errors(status_convert_to_csc);

	//        convert via sparseBLAS API to Handle containing internal data for
	//        subsequent Inspector-executor Sparse BLAS operations.
	status_to_csr = mkl_sparse_s_create_csr ( &quantity_matrix , SPARSE_INDEX_BASE_ZERO,
			quantity_rows, quantity_columns, quantity_IA, quantity_IA+1, quantity_JA, quantity_csr_values );
	check_errors(status_to_csr);
	/** ---------------------------------------------------------------------------
	 ** Populate Shipdate Matrix
	 ** -------------------------------------------------------------------------*/
	//read shipdate
	tbl_read(
			table_file , 11, &shipdate_nnz, &shipdate_rows, &shipdate_columns ,
			&shipdate_csr_values, &shipdate_JA, &shipdate_IA
		);

	// Memory Allocation
	shipdate_csc_values = (float*) mkl_malloc (( shipdate_nnz * sizeof(float) ), MEM_LINE_SIZE );
	shipdate_JA1_csc = (MKL_INT*) mkl_malloc (( shipdate_nnz * sizeof(MKL_INT) ), MEM_LINE_SIZE );
	shipdate_IA1_csc = (MKL_INT*) mkl_malloc (((shipdate_nnz+1) * sizeof(MKL_INT)), MEM_LINE_SIZE );

	// Convert from CSR to CSC
	mkl_scsrcsc(job_csr_csc, &shipdate_nnz, shipdate_csr_values, shipdate_JA, shipdate_IA, shipdate_csc_values, shipdate_JA1_csc, shipdate_IA1_csc, &status_convert_to_csc);
	printf("conversion of shipdate matrix from CSR to CSC ok?\n\t");
	check_errors(status_convert_to_csc);
	//        convert via sparseBLAS API to Handle containing internal data for
	//        subsequent Inspector-executor Sparse BLAS operations.
	status_to_csr = mkl_sparse_s_create_csr (
			&shipdate_matrix , SPARSE_INDEX_BASE_ZERO,      shipdate_rows, shipdate_columns,
			shipdate_IA, shipdate_IA+1, shipdate_JA, shipdate_csr_values
			);

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
	bang_vector = (float*) malloc ( ((quantity_columns+1) * sizeof(float)));
	for (int pos =0; pos < line_status_columns ; pos++){
		bang_vector[pos] = 1.0;
	}

	aggregation_vector = (float*) mkl_malloc ( ((quantity_columns+1) * sizeof(float)), MEM_LINE_SIZE );
	intermediate_vector = (float*) mkl_malloc ( ((quantity_columns+1) * sizeof(float)), MEM_LINE_SIZE );
	final_vector = (float*) mkl_malloc ( ((quantity_columns+1) * sizeof(float)), MEM_LINE_SIZE );

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

	////////////////////////
	// START TIME MEASUREMENT
	////////////////////////
printf("** START TIME MEASUREMENT\n");	
GET_TIME(global_time_start);

	csc_to_csc_mx_selection_and(
			shipdate_csc_values, shipdate_JA1_csc, shipdate_IA1_csc,
			shipdate_nnz, shipdate_rows, shipdate_columns,
			GREATER_EQ , "1998-08-28", LESS_EQ , "1998-12-01",
			&selection_csc_values, &selection_JA1_csc, &selection_IA1_csc,
			&selection_nnz, &selection_rows, &selection_columns
			);

	GET_TIME(global_time_selection);
printf("projection  lw\n");
	// compute projection_lw = line_status krao selection
	csc_csc_krao(
			line_status_csc_values, line_status_JA1_csc, line_status_IA1_csc,
			line_status_nnz, line_status_rows, line_status_columns,
			selection_csc_values, selection_JA1_csc, selection_IA1_csc ,
			selection_nnz, selection_rows, selection_columns,
			&projection_lw_csc_values, &projection_lw_JA1_csc, &projection_lw_IA1_csc,
			&projection_lw_nnz, &projection_lw_rows, &projection_lw_columns
		    );

	GET_TIME(global_time_projection_lw);


printf("projection  rq\n");
	// compute projection_rq = return_flag krao quantity
	csc_csc_krao(
			return_flag_csc_values, return_flag_JA1_csc, return_flag_IA1_csc,
			return_flag_nnz, return_flag_rows, return_flag_columns,
			quantity_csc_values, quantity_JA1_csc, quantity_IA1_csc ,
			quantity_nnz, quantity_rows, quantity_columns,
			&projection_rq_csc_values, &projection_rq_JA1_csc, &projection_rq_IA1_csc,
			&projection_rq_nnz, &projection_rq_rows, &projection_rq_columns
		    );

	GET_TIME(global_time_projection_rq);

printf("projection  final\n");
	//compute projection_final = projection_rq krao projection_lw
	csc_csr_krao(
			projection_rq_csc_values, projection_rq_JA1_csc, projection_rq_IA1_csc,
			projection_rq_nnz, projection_rq_rows, projection_rq_columns,

			projection_lw_csc_values, projection_lw_JA1_csc, projection_lw_IA1_csc,
			projection_lw_nnz, projection_lw_rows, projection_lw_columns,

			&projection_final_csr_values, &projection_final_JA, &projection_final_IA,
			&projection_final_nnz, &projection_final_rows, &projection_final_columns
		    );

	//        convert via sparseBLAS API to Handle containing internal data for
	//        subsequent Inspector-executor Sparse BLAS operations.
	status_to_csr = mkl_sparse_s_create_csr (
			&projection_final_matrix , SPARSE_INDEX_BASE_ZERO, projection_final_rows, projection_final_columns,
			projection_final_IA, projection_final_JA+1, projection_final_JA, projection_final_csr_values
			);
	GET_TIME(global_time_projection_final);

printf("compute final result\n");

	// compute final_result = projection_final * bang
	//
	final_result = mkl_sparse_s_mv (
			SPARSE_OPERATION_NON_TRANSPOSE, 1.0, projection_final_matrix , descrA, bang_vector, 0.0,  final_vector
			);

	////////////////////////
	// STOP TIME MEASUREMENT
	////////////////////////
	GET_TIME(global_time_stop);
	printf("** STOP TIME MEASUREMENT\n");
	writeResults( argv[1] );
	for (int pos =0; pos < projection_final_columns ; pos++){
		if ( final_vector[pos] != 0) {
			printf("%f \n", final_vector[pos]);
		}
	}
	return 0;
}

