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
  int* matrix_values;
  long* matrix_rows;
  long* column_index;
  int current_values_size = ARRAY_SIZE;
  matrix_values = (int*) malloc (current_values_size * sizeof(int));
  matrix_rows = (long*) malloc (current_values_size * sizeof(long));
  column_index = (long*) malloc (current_values_size * sizeof(long));

  FILE* stream = fopen(argv[1], "r");
  char line[1024];
  int column = atoi (argv[2]);
  int value_num = 0;
  while (fgets(line, 1024, stream))
  {
    char *field = (char*) malloc( 128 * sizeof(char) );
    char* tmp = strdup(line);
    field = getfield(tmp, column, field);
    //GQuark* quark;
    GQuark quark = g_quark_from_string (field);

    if ( value_num > current_values_size ){
      matrix_values = (int*) realloc(matrix_values, current_values_size * GROWTH_FACTOR * sizeof(int) );
      matrix_rows =  (long*) realloc(matrix_rows, current_values_size * GROWTH_FACTOR * sizeof(long) );
      column_index =  (long*) realloc(column_index, current_values_size * GROWTH_FACTOR * sizeof(long) );
      current_values_size *= GROWTH_FACTOR;
    }
    matrix_values[value_num]=1;
    //matrix_rows[value_num]=quark;
    //column_index[value_num]=quark;


    printf("%d ,%d \n", value_num, quark);
    // printf("%d - %d \t Field %d would be: {%s}\tQUARK: {%d}\tRE-CONVERTED: {%s}\n", line_num, quark, column, field, quark, gnu_string );
    free(tmp);
    value_num++;
  }
  printf("total columns %d\n int: %dquark: %d\n", value_num, sizeof(int), sizeof(GQuark*));
  return 0;
}
