
// system
#include <iostream>

// libmesh
#include "libmesh/libmesh_logging.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/dense_vector.h"

// hack to allow access to private/protected members
#define private public
#define protected public

// channel
#include "channel_config.h"
#include "spalart_allmaras.h"

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

  GetPot inputfile;
  inputfile = GetPot( );
  AGNOS::PhysicsChannelFlow<T_S,T_P> flowSolver(
      Communicator(MPI_COMM_NULL),inputfile
      );

  BOOST_REQUIRE( flowSolver._availableSolutions.count("primal") );
  BOOST_REQUIRE( flowSolver._availableSolutions.count("adjoint") );
  BOOST_REQUIRE( flowSolver._availableSolutions.count("qoi") );
  BOOST_REQUIRE( flowSolver._availableSolutions.count("errorEstimate") );
  BOOST_REQUIRE( flowSolver._availableSolutions.count("errorIndicators") );
  
  // check for initial parameters set to NULL values
  BOOST_REQUIRE( &flowSolver.getSystem() != NULL );
  BOOST_REQUIRE( &flowSolver.getEquationSystems() != NULL );
  BOOST_REQUIRE( &flowSolver.getMesh() != NULL );
  BOOST_REQUIRE( &flowSolver.getMeshRefinement() != NULL );
  BOOST_REQUIRE( &flowSolver.getEstimator() != NULL );
  BOOST_REQUIRE( &flowSolver.getQois() != NULL );
}

BOOST_AUTO_TEST_CASE( channel_solve )
{
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];

  char* name = new char[27];
  strcpy(name, "dummy");
  
  av[0] = name;
  
  // Initialize libmesh
  LibMeshInit init (ac, av);

  GetPot inputfile;
  inputfile = GetPot( );
  inputfile.set("channel_input","flow.in") ;
  AGNOS::PhysicsChannelFlow<T_S,T_P> flowSolver(
      Communicator(MPI_COMM_NULL),inputfile
      );

  // test setting parameter values
  std::vector<double> params(7,1.);
  T_S parameterValues(params) ;
  flowSolver._setParameterValues(parameterValues);
  ChannelSystem* system = static_cast<ChannelSystem*>(
      flowSolver._system ) ;
  spalartAllmarasTurbulenceModel* turb =
    static_cast<spalartAllmarasTurbulenceModel*>(
        system->_turbModel.get() );

  BOOST_CHECK_CLOSE( turb->_kap, params[0]*0.41 ,1e-9);
  BOOST_CHECK_CLOSE( turb->_cb1 , params[1]*0.1355, 1e-9);
  BOOST_CHECK_CLOSE( turb->_sig  , params[2]*2.0/3.0, 1e-9);
  BOOST_CHECK_CLOSE( turb->_cb2    , params[3]*0.622, 1e-9);
  BOOST_CHECK_CLOSE( turb->_cv1 , params[4]*7.1, 1e-9);
  BOOST_CHECK_CLOSE( turb->_cw2 , params[5]*0.3, 1e-9);
  BOOST_CHECK_CLOSE( turb->_cw3 , params[6]*2.0, 1e-9);

  // make sure nothing catastrophic happens if we try to solve
  flowSolver._solve() ;

}
