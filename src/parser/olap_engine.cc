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

#include "olap_engine.hxx"

int OLAP::OLAP_Engine::get_row_from_string(std::string field ){
  return (((int) g_quark_from_string ( field.c_str() )) -1 );
}

