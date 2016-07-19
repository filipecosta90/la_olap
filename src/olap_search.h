/* ---------------------------------------------------------------------------
 **    Filename: olap_search.h
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
 **
 **     Authors: Filipe Oliveira <a57816@alunos.uminho.pt>
 **          and Sérgio Caldas   <a57779@alunos.uminho.pt>
 **
 ** University of Minho, High Performance Computing Dpt. , April 2016
 ** -------------------------------------------------------------------------*/

//Cache-Lines size is (typically) 32 bytes
#define MEM_LINE_SIZE 32
#define ARRAY_SIZE MEM_LINE_SIZE 
#define GROWTH_FACTOR 2
#define MAX_FIELD_SIZE 128
#define MAX_REG_SIZE 1024
#define LESS 1
#define LESS_EQ 2
#define GREATER 3
#define GREATER_EQ 4

#ifndef _olap_h
#define _olap_h

////////////////////////////////////// AUX /////////////////////////////////////
////////////////////////////////////// AUX /////////////////////////////////////
////////////////////////////////////// AUX /////////////////////////////////////

//starts at position 1 
char* getfield( char* line, int num, char* return_string );

void print_csc(
    float* csc_values, int* A_row_ind, int* A_col_ptr, 
    int NNZ, int number_rows, int number_columns 
    );

void print_csc_vector(
    float* csc_values, int* row_ind, 
    int NNZ, int number_rows
    );

void print_csr( 
    float* csr_values, int* JA, int* IA, 
    int NNZ, int number_rows, int number_columns 
    );

void convert_and_write_to_csv (
    char* filename,
    float* csr_values, int* JA, int* IA,
    int NNZ, int number_rows, int number_columns
    );

void read_from_mx (
    char* filename,
    float** A_csr_values, int** A_JA, int** A_IA,
    int* nnz, int* rows, int* columns
    );

void tbl_read(
    char* table_name, int tbl_column,
    int* nnz, int* rows, int* columns,
    float** A_csr_values, int** A_JA, int** A_IA
    );

int tbl_get_number_elements (char* table_name);

void tbl_read_csc (
    char* table_name, int tbl_column, int number_elements,
    int* nnz, int* rows, int* columns,
    float** A_csr_values, int** A_row_ind, int** A_col_ptr
    );

void tbl_read_csc_measure (
    char* table_name, int tbl_column, int number_elements,
    int* nnz, int* rows, int* columns,
    float** A_csr_values, int** A_row_ind, int** A_col_ptr
    );

void col_read_csc (
    char* table_name, int number_elements,
    int* n_nnz, int* n_rows, int* n_cols,
    float** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_csc_values,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_row_ind,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_col_ptr
    );

void col_read_csc_measure (
    char* table_name, int number_elements,
    int* nnz, int* rows, int* columns,
    float** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_csc_values,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_JA,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_IA
    );

void tbl_read_measure(
    char* table_name, int tbl_column,
    int* nnz, int* rows, int* columns,
    float** A_csr_values, int** A_JA, int** A_IA
    );

void tbl_read_filter(
    char* table_name, int tbl_column, int opp_code, char* comparation_key,
    int* nnz, int* rows, int* columns, 
    float** A_csr_values, int** A_JA, int** A_IA,
    int **quark_start_end, int* quark_global_pos
    );

void csr_csr_square_reshape (
    float** A_csr_values, int** A_JA, int** A_IA,
    int *A_n_nnz, int *A_rows, int *A_columns,
    int reshape_square
    );

void check_errors( sparse_status_t stat );

void csc_tbl_write(
    char* table_name,
    float* A_csc_values, int* A_JA1, int* A_IA1,
    int A_NNZ, int A_number_rows, int A_number_columns
    );

void csr_tbl_write(
    char* table_name,
    float* A_csr_values, int* A_JA, int* A_IA,
    int A_NNZ, int A_number_rows, int A_number_columns
    );

void csr_measure_tbl_write(
    char* table_name,
    float* A_csr_values, int* A_JA, int* A_IA,
    int A_NNZ, int A_number_rows, int A_number_columns
    );

void csr_vector_write(
    char* vector_name,
    float* Vector_csr_values, int Vector_NNZ
    );

void csr_measure_vector_write(
    char* vector_name,
    float* Vector_csr_values, int Vector_NNZ
    );

void csc_to_csr_mx_selection_and(
    float* A_csc_values, int* A_JA1, int* A_IA1,
    int A_NNZ, int A_number_rows, int A_number_columns,
    int opp_code, char* comparation_key, int opp_code2, char* comparation_key2,
    float** C_csr_values, int** C_JA, int** C_IA,
    int* C_NNZ, int* C_number_rows, int* C_number_columns
    );

void csc_to_csc_mx_selection_and(
    float* A_csc_values, int* A_row_ind, int* A_col_ptr,
    int A_NNZ, int A_number_rows, int A_number_columns,
    int opp_code, char* comparation_key, int opp_code2, char* comparation_key2,
    float** C_csc_values, int** C_row_ind, int** C_col_ptr,
    int* C_NNZ, int* C_number_rows, int* C_number_columns
    );

void csc_csc_mx_selection_and(
    float* __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_csc_values,
    int* __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_row_ind,
    int* __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_col_ptr,
    int A_NNZ, int A_number_rows, int A_number_columns,
    int opp_code, char *restrict comparation_key, int opp_code2, char*restrict comparation_key2,
    float** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_csc_values,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_row_ind,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_col_ptr,
    int* C_n_nnz, int* C_n_rows, int* C_n_cols
    );

void csr_mx_selection_or(
    float* A_csr_values, int* A_JA, int* A_IA,
    int A_NNZ, int A_number_rows, int A_number_columns,
    int opp_code, char* comparation_key, int opp_code2, char* comparation_key2,
    float** C_csr_values, int** C_JA, int** C_IA,
    int* C_NNZ, int* C_number_rows, int* C_number_columns,
    int **quark_start_end, int* quark_global_pos
    );

void csr_mx_selection(
    float* A_csr_values, int* A_JA, int* A_IA,
    int A_NNZ, int A_number_rows, int A_number_columns,
    int opp_code, char* comparation_key,
    float** C_csr_values, int** C_JA, int** C_IA,
    int* C_NNZ, int* C_number_rows, int* C_number_columns,
    int **quark_start_end, int* quark_global_pos
    );

///////////////////////////////////// OPS /////////////////////////////////////
///////////////////////////////////// OPS /////////////////////////////////////
///////////////////////////////////// OPS /////////////////////////////////////

/////////////////////////////////
//
//   COMPUTE HADAMARD
//
/////////////////////////////////
void csr_hadamard(
    float *restrict A_csr_values, int *restrict A_JA, int *restrict A_IA,
    int A_NNZ, int A_number_rows, int A_number_columns,
    float *restrict B_csr_values, int *restrict B_JA, int *restrict B_IA,
    int B_NNZ, int B_number_rows, int B_number_columns,
    float **restrict C_csr_values, int **restrict C_JA, int **restrict C_IA,
    int* C_NNZ, int* C_number_rows, int* C_number_columns
    );

/////////////////////////////////
//
//   COMPUTE KHATRI-RAO
//
/////////////////////////////////
void csr_csr_krao(
    float *restrict A_csr_values, int *restrict A_JA, int *restrict A_IA,
    int A_NNZ, int A_number_rows, int A_number_columns,
    float *restrict B_csr_values, int *restrict B_JA, int *restrict B_IA ,
    int B_NNZ, int B_number_rows, int B_number_columns,
    float **restrict C_csr_values, int **restrict C_JA, int **restrict C_IA,
    int* C_NNZ, int* C_number_rows, int* C_number_columns
    );

void csc_to_csr_and_csc_krao(
    float *restrict A_csc_values, int *restrict A_row_ind, int *restrict A_col_ptr,
    int A_NNZ, int A_number_rows, int A_number_columns,
    float *restrict B_csc_values, int *restrict B_JA1, int *restrict B_IA1 ,
    int B_NNZ, int B_number_rows, int B_number_columns,
    float **restrict C_csr_values, int **restrict C_JA, int **restrict C_IA,
    float **restrict C_csc_values, int **restrict C_JA1, int **restrict C_IA1,
    int* C_NNZ, int* C_number_rows, int* C_number_columns
    );

void csc_csc_krao(
    float *restrict A_csc_values, int *restrict A_row_ind, int *restrict A_col_ptr,
    int A_n_nnz, int A_n_rows, int A_n_cols,
    float *restrict B_csc_values, int *restrict B_row_ind, int *restrict B_col_ptr,
    int B_n_nnz, int B_n_rows, int B_n_cols,
    float **C_csc_values, int **C_row_ind, int **C_col_ptr,
    int* C_n_nnz, int* C_n_rows, int *C_n_cols
    );

/////////////////////////////////
//
//   COMPUTE KRONECKER PRODUCT
//
/////////////////////////////////

void csr_kron(
    float* A_csr_values, int* A_JA, int* A_IA,
    int A_NNZ, int A_number_rows, int A_number_columns,
    float* B_csr_values, int* B_JA, int* B_IA ,
    int B_NNZ, int B_number_rows, int B_number_columns,
    float** C_csr_values, int** C_JA, int** C_IA,
    int* C_NNZ, int* C_number_rows, int* C_number_columns
    );

void csc_csc_mm(
    float *restrict A_csc_values, int *restrict A_row_ind, int *restrict A_col_ptr,
    int A_n_nnz, int A_n_rows, int A_n_cols,
    float *restrict B_csc_values, int *restrict B_row_ind, int *restrict B_col_ptr,
    int B_n_nnz, int B_n_rows, int B_n_cols,
    float **C_csc_values, int **C_row_ind, int **C_col_ptr,
    int *C_n_nnz, int *C_n_rows, int *C_n_cols
    );

void csc_bang(
    float *restrict A_csc_values, int *restrict A_row_ind, int *restrict A_col_ptr,
    int A_n_nnz, int A_n_rows, int A_n_cols,
    float **C_csc_values, int **C_row_ind,
    int *C_n_nnz, int *C_n_rows
    );

void produce_tuple_from_krao_csc(
    float *restrict C_csc_values, int *restrict C_row_ind, 
    int C_n_nnz, int C_n_rows, 
    int A_n_rows, int B_n_rows    
    );

void write_coo_from_csc(
    char* filename,
    int C_n_nnz, int C_n_rows, int C_n_cols,
    float* C_csc_values, int*  C_row_ind, int*  C_col_ptr
    );
 
#endif

