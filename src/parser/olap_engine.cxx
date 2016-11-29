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

OLAP::OLAP_Engine::OLAP_Engine(){
  cubes_table = g_hash_table_new_full(
      g_str_hash, g_str_equal, //< This is an string hash.
      free, //< Call "free" on the key (made with "malloc").
      free //< Call "free" on the value (made with "strdup").
      );
}

bool OLAP::OLAP_Engine::create_cube( std::string cube_name  ){
  OLAP::OLAP_Cube *cube = new OLAP::OLAP_Cube( cube_name );
  g_hash_table_insert (cubes_table, (void*) cube_name.c_str(), cube);
  return true;
}

OLAP::OLAP_Cube* OLAP::OLAP_Engine::cube_lookup( std::string identifier ){
  OLAP::OLAP_Cube *cube = nullptr;
  cube = (OLAP::OLAP_Cube*) g_hash_table_lookup ( cubes_table,  (void*) identifier.c_str() );
  return cube;

}

int OLAP::OLAP_Engine::get_row_from_string(std::string field ){
  // since quarks start by 1 and we want to start at line 0 and not 1 lets decrement
  return (((int) g_quark_from_string ( field.c_str() )) -1 );
}

