#ifndef __OLAP_ENGINE_HXX__
#define __OLAP_ENGINE_HXX__ 1

#include <glib.h>
#include "olap_cube.hxx"

namespace OLAP{
  class OLAP_Engine {
    public:
      OLAP_Engine(); // constructor
      int get_row_from_string( std::string field );
      bool create_cube( std::string cube_name  );
      OLAP::OLAP_Cube* cube_lookup( std::string identifier ); 

    private:
      GHashTable *cubes_table = NULL;

  };

}/* end namespace OLAP */
#endif /* END __OLAP_ENGINE_HXX__ */
