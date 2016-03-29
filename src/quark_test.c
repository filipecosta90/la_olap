#include <stdio.h>
#include <glib.h>
#include <stdlib.h>
#include <string.h>

char* getfield(char* line, int num, char* return_string ){
  return_string = strtok(line, ",");
  int pos = 1;
  for ( ; pos <= num; pos++ ){
    return_string = strtok(NULL, ",\n");
  }
  return return_string;
}

int main( int argc, char* argv[]){
  GList* list = NULL;
  FILE* stream = fopen(argv[1], "r");

  char line[1024];

  int column = atoi (argv[2]);
  int line_num = 0;
  while (fgets(line, 1024, stream))
  {
    char *field = (char*) malloc( 128 * sizeof(char) );
    char* tmp = strdup(line);
    field = getfield(tmp, column, field);
    GQuark* quark;
    quark = g_quark_from_string (field);
    quark = g_quark_from_string (field);
    gchar* gnu_string;
    gnu_string = g_quark_to_string(quark);
    printf("%d - %d \t Field %d would be: {%s}\tQUARK: {%d}\tRE-CONVERTED: {%s}\n", line_num, quark, column, field, quark, gnu_string );
    list = g_list_append(list, quark ); 
    free(tmp);
    line_num++;
  }
  printf("The first item is '%d' %s \n", g_list_first(list)->data, g_quark_to_string( g_list_first(list)->data ) );
  return 0;
}
