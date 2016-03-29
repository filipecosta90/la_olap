#include <stdio.h>
#include <glib.h>
#include <stdlib.h>
#include <string.h>
#include "mkl_spblas.h"

//Cache-Lines size is (typically) 64 bytes
#define LINE_SIZE 64
#define ARRAY_SIZE LINE_SIZE / sizeof (int)
#define GROWTH_FACTOR 2

char* getfield(char* line, int num, char* return_string ){
  return_string = strtok(line, ",");
  int pos = 1;
  for ( ; pos <= num; pos++ ){
    return_string = strtok(NULL, ",\n");
  }
  return return_string;
}

int main( int argc, char* argv[]){
  //define sparse-matrix M
  float* matrix_values;
  int* matrix_rows;
  int* column_index;
  int current_values_size = ARRAY_SIZE;
  matrix_values = (float*) malloc (current_values_size * sizeof(float));
  matrix_rows = (int*) malloc (current_values_size * sizeof(int));
  column_index = (int*) malloc (current_values_size * sizeof(int));



  FILE* stream = fopen(argv[1], "r");
  char line[1024];
  int column = atoi (argv[2]);

  int key_col = atoi (argv[2]);.
    int column = atoi (argv[3]);.

    int column_num = 1;
  int number_rows = 1;
  for( ; (fgets(line, 1024, stream) ) ; )
  {
    char *key = (char*) malloc( 128 * sizeof(char) );
    char *field = (char*) malloc( 128 * sizeof(char) );
    char* tmp_key = strdup(line);
    char* tmp_field = strdup(line);
    key = getfield(tmp_key, key_col, key);
    field = getfield(tmp_field, column, field);
    GQuark* quark_key;
    quark_key = g_quark_from_string (key);
    GQuark* quark_field;
    quark_field = g_quark_from_string (field);
    if ( column_num > current_values_size ){
      current_values_size *= GROWTH_FACTOR;
      matrix_values = realloc(matrix_values, current_values_size * GROWTH_FACTOR * sizeof(float) );
      matrix_rows =  realloc(matrix_rows, current_values_size * GROWTH_FACTOR * sizeof(int) );
      column_index =  realloc(column_index, current_values_size * GROWTH_FACTOR * sizeof(int) );
    }

    matrix_values[column_num]=1.0;
    matrix_rows[column_num]= (int) quark;
    column_index[column_num]= (int) quark;

    if( quark > number_rows ) {
      number_rows = quark;
    }

    printf("%d ,%d\n", (int) quark , column_num);
    free(tmp);
    column_num++; 
  }

  sparse_matrix_t A;
  int stat;
  stat = mkl_sparse_s_create_csc(&A,SPARSE_INDEX_BASE_ZERO, number_rows, column_num, column_index, column_index +1, column_index, matrix_values );
  if ( stat == SPARSE_STATUS_SUCCESS){
    //  printf( "SPARSE_STATUS_SUCCESS\n");
  }
  //printf("matrix %d >< %d \n", number_rows, column_num); 
  return 0;
}
