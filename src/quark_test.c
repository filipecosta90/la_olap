#include <stdio.h>
#include <glib.h>
#include <stdlib.h>
#include <string.h>

const char* getfield(char* line, int num){
const char* tok;
  for (tok = strtok(line, ",");
      tok && *tok;
      tok = strtok(NULL, ",\n"))
  {
    if (!--num)
      return tok;
  }
  return NULL;
}

int main( int argc, char* argv[]){
  GList* list = NULL;

  FILE* stream = fopen(argv[1], "r");

  char line[1024];
  int column = atoi (argv[2]); 
  while (fgets(line, 1024, stream))
  {
    char* tmp = strdup(line);
    list = g_list_append(list, getfield(tmp,column)); 
    printf("Field 3 would be %s\n", getfield(tmp, column));
    // NOTE strtok clobbers tmp
    free(tmp);
  }

  printf("The first item is '%s'\n", g_list_first(list)->data);
  return 0;
}
