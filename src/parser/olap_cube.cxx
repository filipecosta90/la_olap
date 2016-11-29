#include <cctype>
#include <fstream>
#include <cassert>
#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <string>

//GLIB
#include <glib.h>

#include "olap_cube.hxx"

int OLAP::OLAP_Cube::get_row_from_string(std::string field ){
  // since quarks start by 1 and we want to start at line 0 and not 1 lets decrement
  return (((int) g_quark_from_string ( field.c_str() )) -1 );
}

