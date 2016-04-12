#include <stdio.h>
#include <glib.h>
#include <stdlib.h>
#include <string.h>

char* getfield(char* line, int num, char* return_string ){
  return_string = strtok(line, ",");
  for ( int pos = 1; pos <= num; pos++ ){
    return_string = strtok(NULL, ",\n");
  }
  return return_string;
}

int main( int argc, char* argv[]){
  GList* list = NULL;
  FILE* stream = fopen(argv[1], "r");

  char line[1024];

  int key_col = atoi (argv[2]); 
  int column = atoi (argv[3]); 
  while (fgets(line, 1024, stream))
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
    printf("%d,%d\n", quark_key, quark_field );
  }
  return 0;
}
