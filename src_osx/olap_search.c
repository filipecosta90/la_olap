/* ---------------------------------------------------------------------------
 **    Filename: olap_search.c
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

#ifndef _olap_c
#define _olap_c

#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <glib.h>
#include <string.h>
#include "olap_search.h"
#include "omp.h"


////////////////////////////////////// AUX /////////////////////////////////////
////////////////////////////////////// AUX /////////////////////////////////////
////////////////////////////////////// AUX /////////////////////////////////////

//starts at position 1

char* getfield( char* line, int num, char* return_string ){
  return_string = strtok(line, "|\n");
  int pos = 1;
  for ( ; pos < num; pos++ ){
    return_string = strtok(NULL, "|\n");
  }
  return return_string;
}

void print_csc(
    float* csc_values, int* row_ind, int* col_ptr,
    int NNZ, int number_rows, int number_columns
    ){
  printf("N NONZ: %d\t", NNZ);
  printf("N ROWS: %d\t", number_rows);
  printf("N COLS: %d\n", number_columns);
  if (number_columns <= 100){
    printf("CSC VALUES(%llu): [\t", sizeof(csc_values));
    for (int pos = 0; pos < NNZ; pos++){
      printf("%f, ", csc_values[pos]);
    }
    printf("] \n row_ind:\t");
    for (int pos = 0; pos < NNZ; pos++){
      printf("%d, ", row_ind[pos]);
    }
    printf("] \ncol_ptr:\t");
    for (int pos = 0; pos <= number_columns; pos++){
      printf("%d, ", col_ptr[pos]);
    }
    printf("]\n");
  }
}

void print_csc_vector(
    float* csc_values, int* row_ind, 
    int NNZ, int number_rows
    ){
  printf("N NONZ: %d\t", NNZ);
  printf("N ROWS: %d\n", number_rows);
  if (number_rows <= 100){
    printf("CSC values: [\t");
    for (int pos = 0; pos < NNZ; pos++){
      printf("%f, ", csc_values[pos]);
    }
    printf("]\nrow_ind: [\t");
    for (int pos = 0; pos < NNZ; pos++){
      printf("%d, ", row_ind[pos]);
    }
    printf("]\n");
  }
}


void print_csr(
    float* csr_values, int* JA, int* IA,
    int NNZ, int number_rows, int number_columns
    ){

  printf("N NONZ: %d\t", NNZ);
  printf("N ROWS: %d\t", number_rows);
  printf("N COLS: %d\n", number_columns);
  if (number_columns <= 100){
    printf("CSR_VALUES(%llu) = [\t", sizeof(csr_values));
    for (int pos = 0; pos < NNZ; pos++){
      printf("%f, ", csr_values[pos]);
    }
    printf("] \nJA = [\t");

    for (int pos = 0; pos < NNZ; pos++){
      printf("%d, ", JA[pos]);
    }

    printf("] \nIA = [\t");
    for (int pos = 0; pos <= number_rows; pos++){
      printf("%d, ", IA[pos]);
    }
    printf("] \n");
  }
}

int tbl_get_number_elements (char* table_name){

  FILE* stream = fopen(table_name, "r");
  if (stream == NULL){
    exit(EXIT_FAILURE);
  }

  int element_number = 0;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;

  for( element_number = 0 ; ( read = getline(&line, &len, stream) ) > 0  ; ++element_number ){
  }
  free(line);
  fclose(stream);

  return element_number;
}

void tbl_read_csc (
    char* table_name, int tbl_column, int number_elements,
    int* n_nnz, int* n_rows, int* n_cols,
    float** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_csc_values,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_row_ind,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_col_ptr
    ){
#ifdef D_VERBOSE
  printf("going to read column %d\n", tbl_column);
#endif

  //define CSC auxiliar sparse-matrix
  // values will contain only ones
  __declspec(align(MEM_LINE_SIZE)) float* aux_csc_values;
  aux_csc_values = (float*) _mm_malloc ((number_elements) * sizeof(float) , MEM_LINE_SIZE );

  // JA  points to column starts in A
  __declspec(align(MEM_LINE_SIZE)) int* aux_csc_row_ind;
  aux_csc_row_ind = (int*) _mm_malloc ((number_elements) * sizeof(int) , MEM_LINE_SIZE );

  // IA splits the array A into rows
  __declspec(align(MEM_LINE_SIZE))  int* aux_csc_col_ptr;
  aux_csc_col_ptr = (int*) _mm_malloc ((number_elements+1) * sizeof(int) , MEM_LINE_SIZE );

#ifdef D_DEBUGGING
  assert(aux_csc_values != NULL);
  assert(aux_csc_row_ind != NULL);
  assert(aux_csc_col_ptr != NULL);
#endif

  int fd;
  fd = open(table_name, O_RDONLY | O_NONBLOCK );

  FILE* stream = fdopen(fd, "r");
  if (stream == NULL){
    exit(EXIT_FAILURE);
  }
  int number_rows = - 1;

  float value;

  int quark_field;
  int current_major_row;
  current_major_row = 0;
  int row_of_element;

  char *field=NULL; 
  char line[MAX_REG_SIZE];
  char *ret;
  size_t len = 0;
  printf("big cycle\n");
  for(int element_number = 0 ; element_number < number_elements ; ++element_number ){
    if ( element_number == 1000000 ){
      printf("1 mil\n");
    }      
    if ( element_number == 10000000 ){
      printf("10 mil\n");
    }
    if ( element_number == 15000000 ){
      printf("15 mil\n");
    }  
    if ( element_number == 20000000 ){
      printf("20 mil\n");
    } 
    if ( element_number == 50000000 ){
      printf("50 mil\n");
    } 
    if ( element_number == 100000000 ){
      printf("100 mil\n");
    } 
    if ( element_number == 125000000 ){
      printf("125 mil\n");
    } 

    if ( element_number == 150000000 ){
      printf("150 mil\n");
    }


    if ( element_number == 175000000 ){
      printf("175 mil\n");
    } 

    ret = fgets(line, MAX_REG_SIZE, stream);
    if (ret == NULL){
      exit(EXIT_FAILURE);
    }
    field = getfield(line, tbl_column, field);

#ifdef D_DEBUGGING
    assert(field!=NULL);
#endif
    quark_field = (int) g_quark_from_string (field);
    // since quarks start by 1 and we want to start at line 0 and not 1 lets decrement
    row_of_element = quark_field -1;
    // for calculating the number of rows
    if (current_major_row < row_of_element ){
      current_major_row = row_of_element ;
    }
    aux_csc_col_ptr[element_number] = element_number;
    aux_csc_row_ind[element_number] = row_of_element;
    aux_csc_values[element_number] = 1.0f;
  }
  free(line);
  fclose(stream);


  aux_csc_col_ptr[number_elements] = number_elements;

  // will contain only ones
  *A_csc_values = aux_csc_values;
  // JA  points to column starts in A
  *A_col_ptr = aux_csc_col_ptr;
  // IA splits the array A into rows
  *A_row_ind = aux_csc_row_ind;

  *n_rows = (current_major_row + 1);
  *n_cols = number_elements;
  *n_nnz = number_elements;

#ifdef D_VERBOSE
  print_csc(
      *A_csc_values, *A_row_ind, *A_col_ptr,
      *n_nnz, *n_rows, *n_cols
      );
  printf("readed matrix %d >< %d : NNZ %d\n", *n_rows, *n_cols, *n_nnz);
#endif

}


void tbl_read_csc_measure (
    char* table_name, int tbl_column, int number_elements,
    int* nnz, int* rows, int* columns,
    float** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_csc_values,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_JA,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_IA
    ){
#ifdef D_VERBOSE
  printf("going to read column %d\n", tbl_column);
#endif

  //define CSC auxiliar sparse-matrix
  // values will contain only ones
  __declspec(align(MEM_LINE_SIZE)) float* aux_csc_values;
  aux_csc_values = (float*) _mm_malloc (number_elements * sizeof(float) , MEM_LINE_SIZE );

  // JA  points to column starts in A
  __declspec(align(MEM_LINE_SIZE)) int* aux_csc_ja;
  aux_csc_ja = (int*) _mm_malloc ( number_elements * sizeof(int) , MEM_LINE_SIZE );

  // IA splits the array A into rows
  __declspec(align(MEM_LINE_SIZE)) int* aux_csc_ia;
  aux_csc_ia = (int*) _mm_malloc ((number_elements+1) * sizeof(int) , MEM_LINE_SIZE );

#ifdef D_DEBUGGING
  assert(aux_csc_values != NULL);
  assert(aux_csc_ja != NULL);
  assert(aux_csc_ia != NULL);
#endif


  int fd;
  fd = open(table_name, O_RDONLY | O_NONBLOCK );

  FILE* stream = fdopen(fd, "r");
  if (stream == NULL){
    exit(EXIT_FAILURE);
  }

  int element_number = 0;

  float value;

  char *field=NULL;
  char line[MAX_REG_SIZE];
  char *ret;
  size_t len = 0;

  for(int element_number = 0 ; element_number < number_elements ; ++element_number ){
    ret = fgets(line, MAX_REG_SIZE, stream);
    if (ret == NULL){
      exit(EXIT_FAILURE);
    }

    field = getfield(line, tbl_column, field);
    value = atof(field);
    aux_csc_ja[element_number] = element_number;
    aux_csc_ia[element_number] = element_number;
    aux_csc_values[element_number] = value;
  }
  free(line);

  fclose(stream);


  aux_csc_ja[number_elements] = number_elements;

  // will contain only ones
  *A_csc_values = aux_csc_values;
  // JA  points to column starts in A
  *A_JA = aux_csc_ja;
  // IA splits the array A into rows
  *A_IA = aux_csc_ia;

  *rows = number_elements;
  *columns = number_elements;
  *nnz = number_elements;

#ifdef D_VERBOSE
  print_csc(
      *A_csc_values, *A_JA, *A_IA,
      *nnz, *rows, *columns
      );
  printf("readed matrix %d %d : NNZ %d\n", *rows, *columns, *nnz);
#endif

}




void col_read_csc (
    char* table_name, int number_elements,
    int* n_nnz, int* n_rows, int* n_cols,
    float** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_csc_values,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_row_ind,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_col_ptr
    ){

  //define CSC auxiliar sparse-matrix
  // values will contain only ones
  __declspec(align(MEM_LINE_SIZE)) float* aux_csc_values;
  aux_csc_values = (float*) _mm_malloc ((number_elements) * sizeof(float) , MEM_LINE_SIZE );

  // JA  points to column starts in A
  __declspec(align(MEM_LINE_SIZE)) int* aux_csc_row_ind;
  aux_csc_row_ind = (int*) _mm_malloc ((number_elements) * sizeof(int) , MEM_LINE_SIZE );

  // IA splits the array A into rows
  __declspec(align(MEM_LINE_SIZE))  int* aux_csc_col_ptr;
  aux_csc_col_ptr = (int*) _mm_malloc ((number_elements+1) * sizeof(int) , MEM_LINE_SIZE );

#ifdef D_DEBUGGING
  assert(aux_csc_values != NULL);
  assert(aux_csc_row_ind != NULL);
  assert(aux_csc_col_ptr != NULL);
#endif

  int fd;
  fd = open(table_name, O_RDONLY | O_NONBLOCK );

  FILE* stream = fdopen(fd, "r");
  if (stream == NULL){
    exit(EXIT_FAILURE);
  }
  int number_rows = - 1;

  float value;

  int quark_field;
  int current_major_row;
  current_major_row = 0;
  int row_of_element;

  char field[MAX_REG_SIZE];
  printf("big cycle\n");

  for(int element_number = 0 ; element_number < number_elements ; ++element_number ){
    if ( element_number == 1000000 ){
      printf("1 mil\n");
    }
    if ( element_number == 10000000 ){
      printf("10 mil\n");
    }
    if ( element_number == 15000000 ){
      printf("15 mil\n");
    }
    if ( element_number == 20000000 ){
      printf("20 mil\n");
    }
    if ( element_number == 50000000 ){
      printf("50 mil\n");
    }
    if ( element_number == 100000000 ){
      printf("100 mil\n");
    }
    if ( element_number == 125000000 ){
      printf("125 mil\n");
    }

    if ( element_number == 150000000 ){
      printf("150 mil\n");
    }

    if ( element_number == 175000000 ){
      printf("175 mil\n");
    }

    fgets(field, MAX_REG_SIZE, stream);

#ifdef D_DEBUGGING
    assert(field!=NULL);
#endif

    quark_field = (int) g_quark_from_string (field);
    // since quarks start by 1 and we want to start at line 0 and not 1 lets decrement
    row_of_element = quark_field -1;
    // for calculating the number of rows
    if (current_major_row < row_of_element ){
      current_major_row = row_of_element ;
    }
    aux_csc_col_ptr[element_number] = element_number;
    aux_csc_row_ind[element_number] = row_of_element;
    aux_csc_values[element_number] = 1.0f;
  }
  fclose(stream);


  aux_csc_col_ptr[number_elements] = number_elements;

  // will contain only ones
  *A_csc_values = aux_csc_values;
  // JA  points to column starts in A
  *A_col_ptr = aux_csc_col_ptr;
  // IA splits the array A into rows
  *A_row_ind = aux_csc_row_ind;

  *n_rows = (current_major_row + 1);
  *n_cols = number_elements;
  *n_nnz = number_elements;

#ifdef D_VERBOSE
  print_csc(
      *A_csc_values, *A_row_ind, *A_col_ptr,
      *n_nnz, *n_rows, *n_cols
      );
  printf("readed matrix %d >< %d : NNZ %d\n", *n_rows, *n_cols, *n_nnz);
#endif

}


void col_read_csc_measure (
    char* table_name, int number_elements,
    int* nnz, int* rows, int* columns,
    float** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_csc_values,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_JA,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_IA
    ){

  //define CSC auxiliar sparse-matrix
  // values will contain only ones
  __declspec(align(MEM_LINE_SIZE)) float* aux_csc_values;
  aux_csc_values = (float*) _mm_malloc (number_elements * sizeof(float) , MEM_LINE_SIZE );

  // JA  points to column starts in A
  __declspec(align(MEM_LINE_SIZE)) int* aux_csc_ja;
  aux_csc_ja = (int*) _mm_malloc ( number_elements * sizeof(int) , MEM_LINE_SIZE );

  // IA splits the array A into rows
  __declspec(align(MEM_LINE_SIZE)) int* aux_csc_ia;
  aux_csc_ia = (int*) _mm_malloc ((number_elements+1) * sizeof(int) , MEM_LINE_SIZE );

#ifdef D_DEBUGGING
  assert(aux_csc_values != NULL);
  assert(aux_csc_ja != NULL);
  assert(aux_csc_ia != NULL);
#endif


  int fd;
  fd = open(table_name, O_RDONLY | O_NONBLOCK );

  FILE* stream = fdopen(fd, "r");
  if (stream == NULL){
    exit(EXIT_FAILURE);
  }

  int element_number = 0;

  float value;

  char field[MAX_REG_SIZE];
  printf("big cycle\n");

  for(int element_number = 0 ; element_number < number_elements ; ++element_number ){
    if ( element_number == 1000000 ){
      printf("1 mil\n");
    }
    if ( element_number == 10000000 ){
      printf("10 mil\n");
    }
    if ( element_number == 15000000 ){
      printf("15 mil\n");
    }
    if ( element_number == 20000000 ){
      printf("20 mil\n");
    }
    if ( element_number == 50000000 ){
      printf("50 mil\n");
    }
    if ( element_number == 100000000 ){
      printf("100 mil\n");
    }
    if ( element_number == 125000000 ){
      printf("125 mil\n");
    }

    if ( element_number == 150000000 ){
      printf("150 mil\n");
    }

    if ( element_number == 175000000 ){
      printf("175 mil\n");
    }

    fgets(field, MAX_REG_SIZE, stream);

    value = atof(field);
    aux_csc_ja[element_number] = element_number;
    aux_csc_ia[element_number] = element_number;
    aux_csc_values[element_number] = value;
  }

  fclose(stream);


  aux_csc_ja[number_elements] = number_elements;

  // will contain only ones
  *A_csc_values = aux_csc_values;
  // JA  points to column starts in A
  *A_JA = aux_csc_ja;
  // IA splits the array A into rows
  *A_IA = aux_csc_ia;

  *rows = number_elements;
  *columns = number_elements;
  *nnz = number_elements;

#ifdef D_VERBOSE
  print_csc(
      *A_csc_values, *A_JA, *A_IA,
      *nnz, *rows, *columns
      );
  printf("readed matrix %d %d : NNZ %d\n", *rows, *columns, *nnz);
#endif

}

void csr_csr_square_reshape (
    float** A_csr_values, int** A_JA, int** A_IA,
    int *A_nnz, int *A_rows, int *A_columns,
    int reshape_square
    ){

  int current_row = (*A_rows);
  int current_column = (*A_columns);
  int current_nnz = (*A_nnz);
  int rows_needed = reshape_square - current_row;
  int columns_needed = reshape_square - current_column;
  int new_nnz = current_nnz + columns_needed;
  int new_rows = reshape_square;
  int new_cols = reshape_square;

  float *new_csr_values;
  int *new_ja;
  int *new_ia; 

#ifdef D_DEBUGGING
  printf("reshaping form %d x %d (%d) to %d x %d (%d)\n", current_row, current_column, current_nnz, new_rows, new_cols, new_nnz );
#endif

  new_csr_values = (float*) malloc ( new_nnz * sizeof(float));
  new_ja = (int*) malloc ( new_nnz * sizeof(int));
  new_ia = (int*) malloc ( (new_rows + 1) * sizeof(int));

  memcpy(new_csr_values, *A_csr_values, current_nnz * sizeof(float));
  memcpy(new_ja, *A_JA, current_nnz * sizeof(int));
  memcpy(new_ia, *A_IA, current_row * sizeof(int));

  for (int at_column = current_column; at_column < new_cols; ++at_column){
    (new_csr_values)[at_column] = 0.0;
    (new_ja)[at_column] = at_column;
  }
  int nnz_it = current_nnz;

  for ( int at_row = current_row ; at_row <= new_rows; ++at_row){
    (new_ia)[at_row] = nnz_it;
    ++nnz_it;
  }

  *A_rows = reshape_square ;
  *A_columns = reshape_square;
  *A_nnz = new_nnz;
  *A_csr_values = new_csr_values;
  *A_JA = new_ja;
  *A_IA = new_ia; 

}



void csc_csc_square_reshape (
    float** A_csc_values, int** A_JA_csc, int** A_IA_csc,
    int *A_nnz, int *A_rows, int *A_columns,
    int reshape_square
    ){

  int current_row = (*A_rows);
  int current_column = (*A_columns);
  int current_nnz = (*A_nnz);


  int rows_needed = reshape_square - current_row;
  int columns_needed = reshape_square - current_column;
  int new_nnz = current_nnz + columns_needed;
  int new_rows = reshape_square;
  int new_cols = reshape_square;
#ifdef D_DEBUGGING
  printf("reshaping form %d x %d (%d) to %d x %d (%d)\n", current_row, current_column, current_nnz, new_rows, new_cols, new_nnz );
#endif
  float* temp_values = (float*) realloc ( (*A_csc_values) , new_nnz * sizeof(float));
  if ( temp_values != NULL ){
    *A_csc_values = temp_values;
  }
  else {
    printf("realloc error in reshape\n");
  }

  int* temp_ja = (int*) realloc ( (*A_JA_csc) , new_nnz * sizeof(int));
  if (temp_ja != NULL ){
    *A_JA_csc = temp_ja;
  }
  else {
    printf("realloc error in reshape\n");
  }

  int* temp_ia = (int*) realloc ( (*A_IA_csc) , (new_nnz + 1) * sizeof(int));
  if (temp_ia != NULL ){
    *A_IA_csc = temp_ia;
  }
  else {
    printf("realloc error in reshape\n");
  }


  for (int at_column = current_column; at_column < new_cols; ++at_column){
    (*A_csc_values)[at_column] = 0.0;
    (*A_JA_csc)[at_column] = at_column;
  }
  int nnz_it = current_nnz;

  for ( int at_row = current_row ; at_row <= new_rows; ++at_row){
    (*A_IA_csc)[at_row] = nnz_it;
    ++nnz_it;
  }

  *A_rows = reshape_square ;
  *A_columns = reshape_square;
  *A_nnz = new_nnz;
}

void csc_to_csc_mx_selection_and(
    float* __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_csc_values,
    int* __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_row_ind,
    int* __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_col_ptr,
    int A_NNZ, int A_number_rows, int A_number_columns,
    int opp_code, char *restrict comparation_key, int opp_code2, char*restrict comparation_key2,
    float** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_csc_values,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_row_ind,
    int** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_col_ptr,
    int* C_n_nnz, int* C_n_rows, int* C_n_cols
    ){

  char* field;

  int non_zero = 0;
  int index;
  int cols;
  int returned_strcmp;
  int returned_strcmp2;
  int at_non_zero = 0;
  int at_row = 0;
  int max_row = 0;

  //define CSC auxiliar sparse-matrix
  // values will contain only ones
  __declspec(align(MEM_LINE_SIZE)) float* aux_csc_values;
  aux_csc_values = (float*) _mm_malloc ((A_NNZ) * sizeof(float) , MEM_LINE_SIZE );

  // JA  points to column starts in A
  __declspec(align(MEM_LINE_SIZE)) int* aux_csc_row_ind;
  aux_csc_row_ind = (int*) _mm_malloc ((A_NNZ) * sizeof(int) , MEM_LINE_SIZE );

  // IA splits the array A into rows
  __declspec(align(MEM_LINE_SIZE)) int* aux_csc_col_ptr;
  aux_csc_col_ptr = (int*) _mm_malloc ((A_number_columns+1) * sizeof(int) , MEM_LINE_SIZE );

#ifdef D_DEBUGGING
  assert(aux_csc_values != NULL);
  assert(aux_csc_row_ind != NULL);
  assert(aux_csc_col_ptr != NULL);
#endif

  __assume_aligned(A_csc_values, MEM_LINE_SIZE);
  __assume_aligned(A_row_ind, MEM_LINE_SIZE);
  __assume_aligned(A_col_ptr, MEM_LINE_SIZE);

  for ( int at_column = 0; at_column < A_number_columns; ++at_column){
    aux_csc_col_ptr[at_column] =  at_non_zero;
    const int iaa = A_row_ind[at_column] +1 ;
    non_zero = 0;
    field = (char*) g_quark_to_string ( iaa );
    returned_strcmp = strcmp( field, comparation_key);
    if (
        ( (opp_code == LESS)  && (returned_strcmp < 0 ))
        ||
        ( (opp_code == LESS_EQ)  && (returned_strcmp <= 0 ))
        ||
        ( (opp_code == GREATER)  && (returned_strcmp > 0 ))
        ||
        ( (opp_code == GREATER_EQ)  && (returned_strcmp >= 0 ))
       ){
      returned_strcmp2 = strcmp( field, comparation_key2);
      if (
          ( (opp_code2 == LESS)  && (returned_strcmp2 < 0) )
          ||
          ( (opp_code2 == LESS_EQ)  && (returned_strcmp2 <= 0) )
          ||
          ( (opp_code2 == GREATER)  && (returned_strcmp2 > 0) )
          ||
          ( (opp_code2 == GREATER_EQ)  && (returned_strcmp2 >= 0) )
         ){
        aux_csc_row_ind[at_non_zero] =  at_column;
        aux_csc_values[at_non_zero] =  A_csc_values[at_column];
        at_non_zero++;
      }
    }
  }
  max_row =  aux_csc_row_ind[at_non_zero];
  aux_csc_col_ptr[A_number_columns] = at_non_zero;
  *C_n_rows = (max_row+1) ;
  *C_n_cols = A_number_columns;
  *C_n_nnz = at_non_zero;
  *C_csc_values = aux_csc_values;
  *C_col_ptr = aux_csc_col_ptr;
  *C_row_ind = aux_csc_row_ind;
}


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
    ){

  char* field;

  int non_zero = 0;
  int index;
  int cols;
  int returned_strcmp;
  int returned_strcmp2;
  int at_non_zero = 0;
  int at_row = 0;
  int max_row = 0;

  //define CSC auxiliar sparse-matrix
  // values will contain only ones
  __declspec(align(MEM_LINE_SIZE)) float* aux_csc_values;
  aux_csc_values = (float*) _mm_malloc ((A_NNZ) * sizeof(float) , MEM_LINE_SIZE );

  // JA  points to column starts in A
  __declspec(align(MEM_LINE_SIZE)) int* aux_csc_row_ind;
  aux_csc_row_ind = (int*) _mm_malloc ((A_NNZ) * sizeof(int) , MEM_LINE_SIZE );

  // IA splits the array A into rows
  __declspec(align(MEM_LINE_SIZE)) int* aux_csc_col_ptr;
  aux_csc_col_ptr = (int*) _mm_malloc ((A_number_columns+1) * sizeof(int) , MEM_LINE_SIZE );


  __declspec(align(MEM_LINE_SIZE)) int* aux_row_bitmap;
  aux_row_bitmap = (int*) _mm_malloc ((A_number_rows+1) * sizeof(int) , MEM_LINE_SIZE );


#ifdef D_DEBUGGING
  assert(aux_csc_values != NULL);
  assert(aux_csc_row_ind != NULL);
  assert(aux_csc_col_ptr != NULL);
  assert(aux_row_bitmap != NULL);
#endif

  __assume_aligned(aux_row_bitmap, MEM_LINE_SIZE);
  for ( int at_row = 0; at_row < A_number_rows; ++at_row){
    aux_csc_values[at_row] =  0;
    non_zero = 0;
    field = (char*) g_quark_to_string ( at_row+1 );
    returned_strcmp = strcmp( field, comparation_key);
    if (
        ( (opp_code == LESS)  && (returned_strcmp < 0 ))
        ||
        ( (opp_code == LESS_EQ)  && (returned_strcmp <= 0 ))
        ||
        ( (opp_code == GREATER)  && (returned_strcmp > 0 ))
        ||
        ( (opp_code == GREATER_EQ)  && (returned_strcmp >= 0 ))
       ){
      returned_strcmp2 = strcmp( field, comparation_key2);
      if (
          ( (opp_code2 == LESS)  && (returned_strcmp2 < 0) )
          ||
          ( (opp_code2 == LESS_EQ)  && (returned_strcmp2 <= 0) )
          ||
          ( (opp_code2 == GREATER)  && (returned_strcmp2 > 0) )
          ||
          ( (opp_code2 == GREATER_EQ)  && (returned_strcmp2 >= 0) )
         ){
        aux_row_bitmap[at_row] =  1;
      }
    }
  }

  __assume_aligned(A_csc_values, MEM_LINE_SIZE);
  __assume_aligned(A_row_ind, MEM_LINE_SIZE);
  __assume_aligned(A_col_ptr, MEM_LINE_SIZE);

  for ( int at_column = 0; at_column < A_number_columns; ++at_column){
    aux_csc_col_ptr[at_column] =  at_non_zero;
    const int iaa = A_row_ind[at_column] +1 ;
    non_zero = 0;
    if ( aux_row_bitmap[iaa] == 1 ){
      aux_csc_row_ind[at_non_zero] =  at_column;
      aux_csc_values[at_non_zero] =  A_csc_values[at_column];
      at_non_zero++;
    }
  }
_mm_free(aux_row_bitmap);
  max_row =  aux_csc_row_ind[at_non_zero];
  aux_csc_col_ptr[A_number_columns] = at_non_zero;
  *C_n_rows = (max_row+1) ;
  *C_n_cols = A_number_columns;
  *C_n_nnz = at_non_zero;
  *C_csc_values = aux_csc_values;
  *C_col_ptr = aux_csc_col_ptr;
  *C_row_ind = aux_csc_row_ind;
}

void csr_mx_selection_or(
    float* A_csr_values, int* A_JA, int* A_IA,
    int A_NNZ, int A_number_rows, int A_number_columns,
    int opp_code, char* comparation_key, int opp_code2, char* comparation_key2,
    float** C_csr_values, int** C_JA, int** C_IA,
    int* C_NNZ, int* C_number_rows, int* C_number_columns,
    int **quark_start_end, int* quark_global_pos
    ){

  *C_csr_values = (float*) malloc ((A_NNZ+1) * sizeof(float));
  *C_JA =  (int*) malloc ((A_NNZ+1)* sizeof(int));
  *C_IA =  (int*) malloc ((A_NNZ+1) * sizeof(int));

  *C_IA[0:A_NNZ] = A_IA[0:A_NNZ];
  *C_JA[0:A_NNZ] = A_JA[0:A_NNZ];
  *C_NNZ = A_NNZ;
  *C_number_rows = A_number_rows;
  *C_number_columns = A_number_columns;

  char* field = (char*) malloc( MAX_FIELD_SIZE * sizeof(char) );

  // read the input file
  for( int element_number = 0 ; element_number < A_NNZ ; ++element_number ){

    int quark_zeroed = 0;

    field = (char*) g_quark_to_string ( element_number );

    int returned_strcmp = strcmp( field , comparation_key );
    int returned_strcmp2 = strcmp( field , comparation_key2 );

    if (
        ( opp_code == LESS  && returned_strcmp >= 0 )
        ||
        ( opp_code == LESS_EQ  && returned_strcmp > 0 )
        ||
        ( opp_code == GREATER  && returned_strcmp <= 0 )
        ||
        ( opp_code == GREATER_EQ  && returned_strcmp < 0 )
       ){
      quark_zeroed = 1;
    }
    else if (
        ( opp_code2 == LESS  && returned_strcmp2 >= 0 )
        ||
        ( opp_code2 == LESS_EQ  && returned_strcmp2 > 0 )
        ||
        ( opp_code2 == GREATER  && returned_strcmp2 <= 0 )
        ||
        ( opp_code2 == GREATER_EQ  && returned_strcmp2 < 0 )
        ){
      quark_zeroed = 1;
    }
    if (quark_zeroed == 1 ){
      *C_csr_values[element_number] = 0;
    }
  }

}


void csr_mx_selection(
    float* A_csr_values, int* A_JA, int* A_IA,
    int A_NNZ, int A_number_rows, int A_number_columns,
    int opp_code, char* comparation_key,
    float** C_csr_values, int** C_JA, int** C_IA,
    int* C_NNZ, int* C_number_rows, int* C_number_columns,
    int **quark_start_end, int* quark_global_pos
    ){


  *C_csr_values = (float*) malloc ((A_NNZ+1) * sizeof(float));
  *C_JA =  (int*) malloc ((A_NNZ+1)* sizeof(int));
  *C_IA =  (int*) malloc ((A_NNZ+1) * sizeof(int));

  *C_IA[0:A_NNZ] = A_IA[0:A_NNZ];
  *C_JA[0:A_NNZ] = A_JA[0:A_NNZ];
  *C_NNZ = A_NNZ;
  *C_number_rows = A_number_rows;
  *C_number_columns = A_number_columns;

  char* field = (char*) malloc( MAX_FIELD_SIZE * sizeof(char) );

  // read the input file
  for( int element_number = 0 ; element_number < A_NNZ ; ++element_number ){

    int quark_zeroed = 0;

    field = (char*) g_quark_to_string ( element_number );

    int returned_strcmp = strcmp( field , comparation_key );

    if (
        ( opp_code == LESS  && returned_strcmp >= 0 )
        ||
        ( opp_code == LESS_EQ  && returned_strcmp > 0 )
        ||
        ( opp_code == GREATER  && returned_strcmp <= 0 )
        ||
        ( opp_code == GREATER_EQ  && returned_strcmp < 0 )
       ){
      quark_zeroed = 1;
    }

    if (quark_zeroed == 1 ){
      *C_csr_values[element_number] = 0;
    }
  }
}

void csc_tbl_write(
    char*  table_name,
    float* A_csc_values, int* A_JA1, int* A_IA1,
    int A_NNZ, int A_number_rows, int A_number_columns
    ){
  char* field = (char*) malloc( MAX_FIELD_SIZE * sizeof(char) );
  FILE* stream = fopen(table_name, "w");
  char line[1024];
  for ( int at_column = 0; at_column < A_number_columns; ++at_column){
    // insert start of column int C_IA1
    int iaa = A_JA1[at_column];
    iaa++;
    field = (char*) g_quark_to_string ( iaa );
    if ( field != NULL  &&  A_csc_values[at_column] > 0 ){
      fprintf(stream, "%s\n", field);
    }
  }
  fclose(stream);
}

void csr_vector_write(
    char* vector_name,
    float* Vector_csr_values, int Vector_NNZ
    ){
#ifdef D_DEBUGGING
  printf("writing vector to file %s\n", vector_name);
#endif

  FILE* stream = fopen(vector_name, "w");
  char* field = (char*) malloc( MAX_FIELD_SIZE * sizeof(char) );

  for ( int at_row = 0; at_row < Vector_NNZ; ++at_row ){
    int iaa = at_row;
    iaa++;
    if (Vector_csr_values[at_row] > 0 ){
      field = (char*) g_quark_to_string ( iaa );
      if ( field != NULL ){
        fprintf(stream, "%s\n", field);
      }
    }
  }
  fclose(stream);
}


void csr_measure_vector_write(
    char* vector_name,
    float* Vector_csr_values, int Vector_NNZ
    ){
#ifdef D_DEBUGGING
  printf("writing vector to file %s\n", vector_name);
#endif

  FILE* stream = fopen(vector_name, "w");

  for ( int at_row = 0; at_row < Vector_NNZ; ++at_row ){
    if (Vector_csr_values[at_row] > 0 ){
      fprintf(stream, "%f\n", Vector_csr_values[at_row]);
    }
  }
  fclose(stream);
}

///////////////////////////////////// OPS /////////////////////////////////////
///////////////////////////////////// OPS /////////////////////////////////////
///////////////////////////////////// OPS /////////////////////////////////////

/////////////////////////////////
//
//   COMPUTE KHATRI-RAO
//
/////////////////////////////////

void csc_csc_krao(
    float *restrict A_csc_values, int *restrict A_row_ind, int *restrict A_col_ptr,
    int A_n_nnz, int A_n_rows, int A_n_cols,
    float *restrict B_csc_values, int *restrict B_row_ind, int *restrict B_col_ptr,
    int B_n_nnz, int B_n_rows, int B_n_cols,
    float **C_csc_values, int **C_row_ind, int **C_col_ptr,
    int *C_n_nnz, int *C_n_rows, int *C_n_cols
    ){

  __assume_aligned(A_csc_values, MEM_LINE_SIZE);
  __assume_aligned(A_row_ind, MEM_LINE_SIZE);
  __assume_aligned(A_col_ptr, MEM_LINE_SIZE);

  __assume_aligned(B_csc_values, MEM_LINE_SIZE);
  __assume_aligned(B_row_ind, MEM_LINE_SIZE);
  __assume_aligned(B_col_ptr, MEM_LINE_SIZE);


  /////////////////////////////////
  //   ALLOCATE MEMORY
  /////////////////////////////////

  __declspec(align(MEM_LINE_SIZE)) float * aux_csc_values;
  __declspec(align(MEM_LINE_SIZE))  int* aux_row_ind;
  __declspec(align(MEM_LINE_SIZE)) int* aux_col_ptr;

  aux_csc_values = (float*) _mm_malloc ( A_n_nnz * sizeof(float) , MEM_LINE_SIZE );
  aux_row_ind = (int*) _mm_malloc ( A_n_nnz  * sizeof(int) , MEM_LINE_SIZE );
  aux_col_ptr = (int*) _mm_malloc ( (A_n_cols+1) * sizeof(int) , MEM_LINE_SIZE );

  /////////////////////////////////
  //   COMPUTE KRAO
  /////////////////////////////////
  const int scalar_B = B_n_rows;
  int current_row = 0;
  int max_row = 0;

#pragma omp parallel
  {
    __assume_aligned(aux_col_ptr, MEM_LINE_SIZE);
    __assume_aligned(A_col_ptr, MEM_LINE_SIZE);
#pragma omp for nowait //reduction(+:aux_col_ptr)
    for ( int at_column = 0 ; at_column < A_n_cols ; ++at_column ){
      aux_col_ptr[at_column] = A_col_ptr[at_column];
    }

    __assume_aligned(A_col_ptr, MEM_LINE_SIZE);
    __assume_aligned(B_col_ptr, MEM_LINE_SIZE);
    __assume_aligned(aux_csc_values, MEM_LINE_SIZE);
    __assume_aligned(A_csc_values, MEM_LINE_SIZE);
    __assume_aligned(B_csc_values, MEM_LINE_SIZE);
#pragma omp for nowait //reduction(+:aux_csc_values)
    for ( int at_column = 0 ; at_column < A_n_cols ; ++at_column ){
      const int a_pos = A_col_ptr[at_column];
      const int b_pos = B_col_ptr[at_column];
      aux_csc_values[a_pos] = A_csc_values[a_pos] * B_csc_values[b_pos];
    }

    __assume_aligned(A_col_ptr, MEM_LINE_SIZE);
    __assume_aligned(B_col_ptr, MEM_LINE_SIZE);
    __assume_aligned(A_row_ind, MEM_LINE_SIZE);
    __assume_aligned(B_row_ind, MEM_LINE_SIZE);
    __assume_aligned(aux_row_ind, MEM_LINE_SIZE);

#pragma omp for nowait //reduction(+:aux_row_ind)
    for ( int at_column = 0 ; at_column < A_n_cols ; ++at_column ){
      const int a_pos = A_col_ptr[at_column];
      const  int b_pos = B_col_ptr[at_column];
      const int current_row = B_row_ind[b_pos] + ( A_row_ind[a_pos] * scalar_B );
      aux_row_ind[a_pos] = current_row;
      if (current_row > max_row){
        max_row = current_row;
      }
    }
  }
  aux_col_ptr[A_n_cols] = A_n_nnz;

  *C_n_rows = (max_row+1);
  *C_n_cols = A_n_cols;
  *C_n_nnz = A_n_nnz;
  *C_csc_values = aux_csc_values;
  *C_row_ind = aux_row_ind;
  *C_col_ptr = aux_col_ptr;

}

/////////////////////////////////
//
//   COMPUTE DOT PRODUCT
//
/////////////////////////////////

void csc_csc_mm(
    float * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_csc_values,
    int * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_row_ind,
    int *__restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_col_ptr,
    int A_n_nnz, int A_n_rows, int A_n_cols,
    float * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) B_csc_values,
    int * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) B_row_ind,
    int * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) B_col_ptr,
    int B_n_nnz, int B_n_rows, int B_n_cols,
    float ** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_csc_values,
    int ** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_row_ind,
    int ** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_col_ptr,
    int *C_n_nnz, int *C_n_rows, int *C_n_cols
    ){

  __assume_aligned(A_csc_values, MEM_LINE_SIZE);
  __assume_aligned(A_row_ind, MEM_LINE_SIZE);
  __assume_aligned(A_col_ptr, MEM_LINE_SIZE);
  __assume_aligned(B_csc_values, MEM_LINE_SIZE);
  __assume_aligned(B_row_ind, MEM_LINE_SIZE);
  __assume_aligned(B_col_ptr, MEM_LINE_SIZE);

  __declspec(align(MEM_LINE_SIZE)) float * aux_csc_values;
  __declspec(align(MEM_LINE_SIZE)) int* aux_row_ind;
  __declspec(align(MEM_LINE_SIZE)) int* aux_col_ptr;
  int b_row = -1;
  int a_row = -1 ;
  int flag_a, flag_b;
  int a_pos, b_pos;
  int max_row = 0;
  int nnz_aux = 0;

  int nnz = A_n_nnz > B_n_nnz ? A_n_nnz : B_n_nnz;

  aux_csc_values = (float*) _mm_malloc ( nnz * sizeof(float) , MEM_LINE_SIZE );
  aux_row_ind = (int*) _mm_malloc ( nnz  * sizeof(int) , MEM_LINE_SIZE);
  aux_col_ptr = (int*) _mm_malloc ( (B_n_cols+1) * sizeof(int) , MEM_LINE_SIZE);

  for ( int at_column_b = 0 ; at_column_b < B_n_cols ; ++at_column_b ){
    aux_col_ptr[at_column_b] = nnz_aux;
    b_pos = B_col_ptr[at_column_b];
    flag_b = B_col_ptr[at_column_b+1] - b_pos;
    if ( flag_b > 0 ) {  
      b_row = B_row_ind[b_pos];
      for ( int at_column_a = 0 ; at_column_a < A_n_cols ; ++at_column_a ){
        a_pos = A_col_ptr[at_column_a]; 
        //check if there is a non zero in this column
        flag_a = A_col_ptr[at_column_a+1] - a_pos;
        if ( ( b_row == at_column_a ) && (flag_a > 0) ){
          a_row = A_row_ind[a_pos];
          aux_row_ind[nnz_aux] = a_row;
          max_row = a_row > max_row ? a_row : max_row;
          aux_csc_values[nnz_aux] += A_csc_values[a_pos] * B_csc_values[b_pos];
          nnz_aux++;
        }
      }
    }
  }
  aux_col_ptr[B_n_cols] = nnz_aux;
  *C_n_rows = (max_row+1);
  *C_n_cols = B_n_cols;
  *C_n_nnz = nnz_aux;
  *C_csc_values = aux_csc_values;
  *C_row_ind = aux_row_ind;
  *C_col_ptr = aux_col_ptr;
}


void csc_csc_bitmap_mm(
    float * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_csc_values,
    int * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_row_ind,
    int *__restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_col_ptr,
    int A_n_nnz, int A_n_rows, int A_n_cols,
    float * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) B_csc_values,
    int * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) B_row_ind,
    int * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) B_col_ptr,
    int B_n_nnz, int B_n_rows, int B_n_cols,
    float ** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_csc_values,
    int ** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_row_ind,
    int ** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_col_ptr,
    int *C_n_nnz, int *C_n_rows, int *C_n_cols
    ){

  __assume_aligned(A_csc_values, MEM_LINE_SIZE);
  __assume_aligned(A_row_ind, MEM_LINE_SIZE);
  __assume_aligned(A_col_ptr, MEM_LINE_SIZE);
  __assume_aligned(B_csc_values, MEM_LINE_SIZE);
  __assume_aligned(B_row_ind, MEM_LINE_SIZE);
  __assume_aligned(B_col_ptr, MEM_LINE_SIZE);

  __declspec(align(MEM_LINE_SIZE)) float * aux_csc_values;
  __declspec(align(MEM_LINE_SIZE)) int* aux_row_ind;
  __declspec(align(MEM_LINE_SIZE)) int* aux_col_ptr;
  int b_row = -1;
  int a_row = -1 ;
  int flag_a, flag_b;
  int a_pos, b_pos;
  int max_row = 0;
  int nnz_aux = 0;

  int nnz = A_n_nnz > B_n_nnz ? A_n_nnz : B_n_nnz;

  aux_csc_values = (float*) _mm_malloc ( nnz * sizeof(float) , MEM_LINE_SIZE );
  aux_row_ind = (int*) _mm_malloc ( nnz  * sizeof(int) , MEM_LINE_SIZE);
  aux_col_ptr = (int*) _mm_malloc ( (B_n_cols+1) * sizeof(int) , MEM_LINE_SIZE);

  for ( int at_column_b = 0 ; at_column_b < B_n_cols ; ++at_column_b ){
    aux_col_ptr[at_column_b] = nnz_aux;
    b_pos = B_col_ptr[at_column_b];
    flag_b = B_col_ptr[at_column_b+1] - b_pos;
    if ( flag_b > 0 ) {
      b_row = B_row_ind[b_pos];
      for ( int at_column_a = 0 ; at_column_a < A_n_cols ; ++at_column_a ){
        a_pos = A_col_ptr[at_column_a];
        flag_a = A_col_ptr[at_column_a+1] - a_pos;
        if ( ( b_row == at_column_a ) && (flag_a > 0) ){
          a_row = A_row_ind[a_pos];
          aux_row_ind[nnz_aux] = a_row;
          max_row = a_row > max_row ? a_row : max_row;
          aux_csc_values[nnz_aux] += 1;
          nnz_aux++;
        }
      }
    }
  }
  aux_col_ptr[B_n_cols] = nnz_aux;
  *C_n_rows = (max_row+1);
  *C_n_cols = B_n_cols;
  *C_n_nnz = nnz_aux;
  *C_csc_values = aux_csc_values;
  *C_row_ind = aux_row_ind;
  *C_col_ptr = aux_col_ptr;
}

void csc_bang(
    float * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_csc_values,
    int * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_row_ind,
    int * __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) A_col_ptr,
    int A_n_nnz, int A_n_rows, int A_n_cols,
    float ** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_csc_values,
    int ** __restrict__  __attribute__((aligned (MEM_LINE_SIZE))) C_row_ind,
    int *C_n_nnz, int *C_n_rows
    ){

  __assume_aligned(A_csc_values, MEM_LINE_SIZE);
  __assume_aligned(A_row_ind, MEM_LINE_SIZE);
  __assume_aligned(A_col_ptr, MEM_LINE_SIZE);

  __declspec(align(MEM_LINE_SIZE)) float * aux_csc_values;
  __declspec(align(MEM_LINE_SIZE)) int* aux_row_ind;
  int a_row = -1 ;
  int max_row = 0;

  int nnz = A_n_rows;
  int nnz_aux = 0 ;
  int flag_found = 0;
  aux_csc_values = (float*) _mm_malloc ( nnz * sizeof(float) , MEM_LINE_SIZE );
  aux_row_ind = (int*) _mm_malloc ( nnz  * sizeof(int) , MEM_LINE_SIZE );

  for ( int at_column_a = 0 ; at_column_a < A_n_cols ; ++at_column_a ){
    const int a_pos = A_col_ptr[at_column_a];
    int flag = A_col_ptr[at_column_a+1] - a_pos;

    if (  flag>0 ){
      a_row = A_row_ind[a_pos];
      flag_found = -1;
      for ( int at_nnz = 0; ( at_nnz <= nnz_aux ) && (flag_found < 0 ); at_nnz++){
        if (aux_row_ind[at_nnz] == a_row){
          flag_found = at_nnz;
        }
      }
      if (flag_found < 0){
        aux_row_ind[nnz_aux] = a_row;
        max_row = a_row > max_row ? a_row : max_row;
        aux_csc_values[nnz_aux] = A_csc_values[a_pos];
        nnz_aux++;
      }
      else {
        aux_csc_values[flag_found] += A_csc_values[a_pos];
      }
    }

  }

  *C_n_rows = (max_row+1);
  *C_n_nnz = nnz_aux;
  *C_csc_values = aux_csc_values;
  *C_row_ind = aux_row_ind;

}

void produce_tuple_from_krao_csc(
    float *restrict C_csc_values, int *restrict C_row_ind, 
    int C_n_nnz, int C_n_rows, 
    int A_n_rows, int B_n_rows    
    ){
  int row_a, row_b, row_c;
  char* field_a;
  char* field_b;

  for (int at_nnz = 0; at_nnz < C_n_nnz; at_nnz++){
    row_c = C_row_ind[at_nnz];
    row_a = row_c / B_n_rows;
    row_b = row_c % B_n_rows;
    row_a++;
    row_b++;
    field_a = (char*) g_quark_to_string ( row_a );
    field_b = (char*) g_quark_to_string ( row_b );
    printf("(%s,%s) %f\n", field_a, field_b,  C_csc_values[at_nnz]);
  }
}

void write_coo_from_csc(
  char* filename,
    int C_n_nnz, int C_n_rows, int C_n_cols,
    float* C_csc_values, int*  C_row_ind, int*  C_col_ptr
    ){
  int row_c;
  char* field_c;
  FILE* stream = fopen(filename, "w");
  if (stream != NULL ){
      fprintf( stream, "row,column\n" ); 
  for (int at_nnz = 0; at_nnz < C_n_nnz; at_nnz++){
    row_c = C_row_ind[at_nnz];
      fprintf( stream, "%d,%d\n", row_c, at_nnz ); 
    }
    fclose(stream);
  }
}

#endif
