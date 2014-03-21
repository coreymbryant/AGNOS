
// system
#include <iostream>

// libmesh
#include "libmesh/libmesh_logging.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/dense_vector.h"

// local
#include "channel_config.h"
#include "discrete_flow.h"

// boost
#define BOOST_TEST_MODULE channelTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


BOOST_AUTO_TEST_CASE( _discrete_flow_interpolate_ )
{
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];

  char* name = new char[27];
  strcpy(name, "_discrete_flow_interpolate_");
  
  av[0] = name;
  
  // Initialize libmesh
  LibMeshInit init (ac, av);

  // turn off performance logging to minimize screen output
  libMesh::perflog.disable_logging();
  
  delete [] name;
  delete [] av;

  // check that libmesh is initialized
  BOOST_CHECK( libMesh::initialized() );
}
