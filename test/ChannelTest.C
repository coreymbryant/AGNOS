
// system
#include <iostream>

// libmesh
#include "libmesh/libmesh_logging.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/dense_vector.h"

// channel
#include "channel_config.h"

// local
#include "PhysicsChannelFlow.h"

// boost
#define BOOST_TEST_MODULE channelTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


BOOST_AUTO_TEST_CASE( channel_constructor )
{
  
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];

  char* name = new char[27];
  strcpy(name, "dummy");
  
  av[0] = name;
  
  // Initialize libmesh
  LibMeshInit init (ac, av);

  {
  GetPot inputfile;
  inputfile = GetPot( );
  inputfile.set("channel_input","flow.in") ;
  AGNOS::PhysicsChannelFlow<T_S,T_P> flowSolver(
      Communicator(MPI_COMM_NULL),inputfile
      );
  }

  BOOST_CHECK( 1 );
}
