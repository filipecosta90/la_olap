/* ---------------------------------------------------------------------------
 **    Filename: olap.h
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
#define ARRAY_SIZE MEM_LINE_SIZE / sizeof (MKL_INT)
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
    float* csc_values, MKL_INT* JA1, MKL_INT* IA1, 
    MKL_INT NNZ, MKL_INT number_rows, MKL_INT number_columns 
    );

void print_csr( 
    float* csr_values, MKL_INT* JA, MKL_INT* IA, 
    MKL_INT NNZ, MKL_INT number_rows, MKL_INT number_columns 
    );

void convert_and_write_to_csv (
    char* filename,
    float* csr_values, MKL_INT* JA, MKL_INT* IA,
    MKL_INT NNZ, MKL_INT number_rows, MKL_INT number_columns
    );

void read_from_csv (
    char* filename,
    float** A_csr_values, MKL_INT** A_JA, MKL_INT** A_IA,
    MKL_INT* nnz, MKL_INT* rows, MKL_INT* columns
    );

void check_errors( sparse_status_t stat );

void tbl_read( 
    char* table_name, MKL_INT tbl_column, 
    MKL_INT* nnz, MKL_INT* rows, MKL_INT* columns, 
    float** A_csr_values, MKL_INT** A_JA, MKL_INT** A_IA
    );

void tbl_read_measure(
    char* table_name, MKL_INT tbl_column,
    MKL_INT* nnz, MKL_INT* rows, MKL_INT* columns,
    float** A_csr_values, MKL_INT** A_JA, MKL_INT** A_IA
    );

void tbl_read_filter( 
    char* table_name, MKL_INT tbl_column, int opp_code, char* comparation_key,
    MKL_INT* nnz, MKL_INT* rows, MKL_INT* columns , 
    float** A_csr_values, MKL_INT** A_JA, MKL_INT** A_IA
    );

void tbl_read_filter_and(
    char* table_name, MKL_INT tbl_column, int opp_code, char* comparation_key, int opp_code2, char* comparation_key2,
    MKL_INT* nnz, MKL_INT* rows, MKL_INT* columns ,
    float** A_csr_values, MKL_INT** A_JA, MKL_INT** A_IA
    );

void tbl_write(
    char* filename, MKL_INT A_nnz, MKL_INT rows, MKL_INT columns,
    float* A_csr_values, MKL_INT* A_JA, MKL_INT* A_IA
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
    float *restrict A_csr_values, MKL_INT *restrict A_JA, MKL_INT *restrict A_IA,
    MKL_INT A_NNZ, MKL_INT A_number_rows, MKL_INT A_number_columns,
    float *restrict B_csr_values, MKL_INT *restrict B_JA, MKL_INT *restrict B_IA,
    MKL_INT B_NNZ, MKL_INT B_number_rows, MKL_INT B_number_columns,
    float **restrict C_csr_values, MKL_INT **restrict C_JA, MKL_INT **restrict C_IA,
    MKL_INT* C_NNZ, MKL_INT* C_number_rows, MKL_INT* C_number_columns
    );

/////////////////////////////////
//
//   COMPUTE KHATRI-RAO
//
/////////////////////////////////
void csr_krao(
    float* A_csr_values, MKL_INT* A_JA, MKL_INT* A_IA, 
    MKL_INT A_NNZ, MKL_INT A_number_rows, MKL_INT A_number_columns,
    float* B_csr_values, MKL_INT* B_JA, MKL_INT* B_IA , 
    MKL_INT B_NNZ, MKL_INT B_number_rows, MKL_INT B_number_columns,
    float** C_csr_values, MKL_INT** C_JA, MKL_INT** C_IA, 
    MKL_INT* C_NNZ, MKL_INT* C_number_rows, MKL_INT* C_number_columns  
    );


/////////////////////////////////
//
//   COMPUTE KRONECKER PRODUCT
//
/////////////////////////////////

void csr_kron(
    float* A_csr_values, MKL_INT* A_JA, MKL_INT* A_IA,
    MKL_INT A_NNZ, MKL_INT A_number_rows, MKL_INT A_number_columns,
    float* B_csr_values, MKL_INT* B_JA, MKL_INT* B_IA ,
    MKL_INT B_NNZ, MKL_INT B_number_rows, MKL_INT B_number_columns,
    float** C_csr_values, MKL_INT** C_JA, MKL_INT** C_IA,
    MKL_INT* C_NNZ, MKL_INT* C_number_rows, MKL_INT* C_number_columns
    );

#endif

