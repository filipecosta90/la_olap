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
 ** University of Minho, High Performance Computing Dpt. , April 2016
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

double global_time_start, global_time_stop, total_time;

void writeResults ( ) {
	total_time = global_time_stop - global_time_start;
	FILE* stream = fopen("timing/timings_no_vec_32.dat", "a+");
	fprintf(stream, "32,%f\n", total_time);
	fclose(stream);
}

int main( int argc, char* argv[]){

	//define CSR sparse-matrix M

	__declspec(align(MEM_LINE_SIZE))  float* returnFlag_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* returnFlag_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* returnFlag_IA;
	MKL_INT returnFlag_rows;
	MKL_INT returnFlag_columns;
	MKL_INT returnFlag_nnz;

	__declspec(align(MEM_LINE_SIZE))  float* lineStatus_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* lineStatus_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* lineStatus_IA;
	MKL_INT lineStatus_rows;
	MKL_INT lineStatus_columns;
	MKL_INT lineStatus_nnz;

	__declspec(align(MEM_LINE_SIZE))  float* quantity_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* quantity_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* quantity_IA;
	MKL_INT quantity_rows;
	MKL_INT quantity_columns;
	MKL_INT quantity_nnz;
	sparse_matrix_t  quantity_matrix;

	__declspec(align(MEM_LINE_SIZE))  float* shipdate_gt_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* shipdate_gt_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* shipdate_gt_IA;
	MKL_INT shipdate_gt_rows;
	MKL_INT shipdate_gt_columns;
	MKL_INT shipdate_gt_nnz;
	sparse_matrix_t  shipdate_gt_matrix;

	__declspec(align(MEM_LINE_SIZE))  float* shipdate_lt_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* shipdate_lt_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* shipdate_lt_IA;
	MKL_INT shipdate_lt_rows;
	MKL_INT shipdate_lt_columns;
	MKL_INT shipdate_lt_nnz;
	sparse_matrix_t  shipdate_lt_matrix;

	__declspec(align(MEM_LINE_SIZE))  float* selection_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* selection_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* selection_IA;
	MKL_INT selection_rows;
	MKL_INT selection_columns;
	MKL_INT selection_nnz;
	sparse_matrix_t  selection_matrix;

	__declspec(align(MEM_LINE_SIZE))  float* projection1_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* projection1_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* projection1_IA;
	MKL_INT projection1_rows;
	MKL_INT projection1_columns;
	MKL_INT projection1_nnz;
	sparse_matrix_t  projection1_matrix;

	__declspec(align(MEM_LINE_SIZE))  float* projection2_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* projection2_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* projection2_IA;
	MKL_INT projection2_rows;
	MKL_INT projection2_columns;
	MKL_INT projection2_nnz;
	sparse_matrix_t  projection2_matrix;

	__declspec(align(MEM_LINE_SIZE)) float* aggregation_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* aggregation_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* aggregation_IA;
	MKL_INT aggregation_rows;
	MKL_INT aggregation_columns;
	MKL_INT aggregation_nnz;

	__declspec(align(MEM_LINE_SIZE)) float* intermediate_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* intermediate_JA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* intermediate_IA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT intermediate_rows;
	MKL_INT intermediate_columns;
	MKL_INT intermediate_nnz;
	sparse_matrix_t  intermediate_matrix;

	__declspec(align(MEM_LINE_SIZE)) float* final_csr_values = NULL;
	__declspec(align(MEM_LINE_SIZE)) MKL_INT* final_JA;
	__declspec(align(MEM_LINE_SIZE)) MKL_INT* final_IA;
	__declspec(align(MEM_LINE_SIZE))  MKL_INT* final_IA_END;
	MKL_INT final_rows;
	MKL_INT final_columns;
	MKL_INT final_nnz;
	sparse_matrix_t  final_matrix;

	__declspec(align(MEM_LINE_SIZE))  float* bang_vector;
	__declspec(align(MEM_LINE_SIZE))  float* aggregation_vector;

	//conversion status from csr arrays into mkl sparse_matrix_t 
	sparse_status_t status_to_csr;

	//read return flag
	read_from_mx("__quark_mx_1/return_flag_1.mx", &returnFlag_csr_values, &returnFlag_JA, &returnFlag_IA, &returnFlag_nnz, &returnFlag_rows, &returnFlag_columns);

	//read line status
	read_from_mx("__quark_mx_1/line_status_1.mx", &lineStatus_csr_values, &lineStatus_JA, &lineStatus_IA, &lineStatus_nnz, &lineStatus_rows, &lineStatus_columns);

	//read quantity
	read_from_mx("__quark_mx_1/quantity_1.mx", &quantity_csr_values, &quantity_JA, &quantity_IA, &quantity_nnz, &quantity_rows, &quantity_columns);

	//        convert via sparseBLAS API to Handle containing internal data for 
	//        subsequent Inspector-executor Sparse BLAS operations.
	status_to_csr = mkl_sparse_s_create_csr ( &quantity_matrix , SPARSE_INDEX_BASE_ZERO, 
			quantity_rows, quantity_columns, quantity_IA, quantity_IA+1, quantity_JA, quantity_csr_values );

	//read shipdate gt
	read_from_mx("__quark_mx_1/shipdate_gt_1.mx", &shipdate_gt_csr_values, &shipdate_gt_JA, &shipdate_gt_IA, &shipdate_gt_nnz, &shipdate_gt_rows, &shipdate_gt_columns);

	//        convert via sparseBLAS API to Handle containing internal data for
	//        subsequent Inspector-executor Sparse BLAS operations.
	status_to_csr = mkl_sparse_s_create_csr ( &shipdate_gt_matrix , SPARSE_INDEX_BASE_ZERO,
			shipdate_gt_rows, shipdate_gt_columns, shipdate_gt_IA, shipdate_gt_IA+1, shipdate_gt_JA, shipdate_gt_csr_values );

	//read shipdate lt
	read_from_mx("__quark_mx_1/shipdate_lt_1.mx", &shipdate_lt_csr_values, &shipdate_lt_JA, &shipdate_lt_IA, &shipdate_lt_nnz, &shipdate_lt_rows, &shipdate_lt_columns);

	//        convert via sparseBLAS API to Handle containing internal data for
	//        subsequent Inspector-executor Sparse BLAS operations.
	status_to_csr = mkl_sparse_s_create_csr ( &shipdate_lt_matrix , SPARSE_INDEX_BASE_ZERO,
			shipdate_lt_rows, shipdate_lt_columns, shipdate_lt_IA, shipdate_lt_IA+1, shipdate_lt_JA, shipdate_lt_csr_values );
	// compute selection = shipdate_gt * shipdate_lt 
	sparse_status_t selection_result;

	////////////////////////
	// START TIME MEASUREMENT
	////////////////////////
	GET_TIME(global_time_start);

	selection_result = mkl_sparse_spmm ( SPARSE_OPERATION_NON_TRANSPOSE,
			shipdate_gt_matrix,
			shipdate_lt_matrix,
			&selection_matrix);

	// compute projection1 = returnFlag krao quantity
	csr_krao(
			returnFlag_csr_values, returnFlag_JA, returnFlag_IA, 
			returnFlag_nnz, returnFlag_rows, returnFlag_columns,
			quantity_csr_values, quantity_JA, quantity_IA ,
			quantity_nnz, quantity_rows,quantity_columns,
			&projection_csr_values, &projection_JA, &projection_IA, 
			&projection_nnz, &projection_rows, &projection_columns  
		);

	status_to_csr = mkl_sparse_s_create_csr ( &projection_matrix , SPARSE_INDEX_BASE_ZERO, projection_rows, projection_columns, projection_IA, projection_IA+1, projection_JA, projection_csr_values );


	// compute projection2 = returnFlag krao quantity
	csr_krao(
			lineStatus_csr_values, lineStatus_JA, lineStatus_IA, 
			lineStatus_nnz, lineStatus_rows, lineStatus_columns,
			selection_csr_values, selection_JA, selection_IA ,
			selection_nnz, selection_rows, selection_columns,
			&projectioni2_csr_values, &projection2_JA, &projection2_IA, 
			&projection2_nnz, &projection2_rows, &projection2_columns  
		);

	//compute aggregation = projection2 * bang

	bang_vector = (float*) mkl_malloc ((quantity_columns * sizeof(float)), MEM_LINE_SIZE );
	aggregation_vector = (float*) mkl_malloc ((quantity_columns * sizeof(float)), MEM_LINE_SIZE );

	sparse_status_t aggregation_result;
	struct matrix_descr descrA;
	descrA.type = SPARSE_MATRIX_TYPE_GENERAL;

	aggregation_result = mkl_sparse_s_mv ( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, projection2_matrix , descrA, bang_vector, 1.0,  aggregation_vector);

	// compute final_result = projection1 * aggregation
	sparse_status_t final_result;

	final_result = mkl_sparse_spmm ( SPARSE_OPERATION_NON_TRANSPOSE, 
			projection1_matrix,
			projection2_matrix,
			&final_matrix);

	////////////////////////
	// STOP TIME MEASUREMENT
	////////////////////////
	GET_TIME(global_time_stop);
	writeResults( );

	return 0;

}
