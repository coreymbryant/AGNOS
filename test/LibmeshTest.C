
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

// local
#include "PhysicsLibmesh.h"

// boost
#define BOOST_TEST_MODULE LibmeshTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_CASE( PhysicsLibmesh_constructor )
{
  
  GetPot inputfile;
  inputfile = GetPot( );
  AGNOS::PhysicsLibmesh<T_S,T_P> physics(
      Communicator(MPI_COMM_NULL),inputfile
      );

  // check for initial parameters set to NULL values
  BOOST_REQUIRE( &physics.getSystem() == NULL );
  BOOST_REQUIRE( &physics.getEquationSystems() == NULL );
  BOOST_REQUIRE( &physics.getMesh() == NULL );
  BOOST_REQUIRE( &physics.getMeshRefinement() == NULL );
  BOOST_REQUIRE( &physics.getEstimator() == NULL );
  BOOST_REQUIRE( &physics.getQois() == NULL );
  BOOST_REQUIRE( physics.comm().get() == MPI_COMM_NULL );
  BOOST_REQUIRE( physics.getAvailableSolutions().count("primal") == 1 );

  // check default values set in constructor
  BOOST_REQUIRE( physics._useUniformRefinement ) ;
  BOOST_REQUIRE( physics._numberHRefinements ==1 ) ;
  BOOST_REQUIRE( physics._numberPRefinements == 0);
  BOOST_REQUIRE( physics._maxRefineSteps == 1);
  BOOST_REQUIRE( !physics._writePrimalViz  );
  BOOST_REQUIRE( !physics._writeAdjointViz );
  BOOST_REQUIRE( !physics._resolveAdjoint );

}

BOOST_AUTO_TEST_CASE( PhysicsLibmesh_initRoutines )
{
  
  GetPot inputfile;
  inputfile = GetPot( );
  AGNOS::PhysicsLibmesh<T_S,T_P> physics(
      Communicator(MPI_COMM_NULL),inputfile
      );
  
  // test buildMeshRefinement 
  physics._mesh = new libMesh::Mesh(physics._communicator);
  physics._buildMeshRefinement();
  BOOST_REQUIRE(physics._meshRefinement-> coarsen_by_parents() );
  BOOST_REQUIRE_CLOSE( 
      physics._meshRefinement->absolute_global_tolerance(),
      1e-6, 1e-16 );
  BOOST_REQUIRE_CLOSE( 
    physics._meshRefinement->refine_fraction(),
    0.7, 1e-16);
  BOOST_REQUIRE_CLOSE( 
    physics._meshRefinement->coarsen_fraction(), 
    0.3, 1e-16 );  
  BOOST_REQUIRE_CLOSE( 
    physics._meshRefinement->coarsen_threshold(),
    1e-5, 1e-16);
  BOOST_REQUIRE( physics._meshRefinement->max_h_level() == 15 );

  // test buildErroEstimator
  physics._qois = new libMesh::QoISet;
  std::vector<unsigned int> qoi_indices;
  qoi_indices.push_back(0);
  physics._qois->add_indices(qoi_indices);
  physics._qois->set_weight(0, 1.0);
  physics._buildErrorEstimator();
  BOOST_REQUIRE(
    physics._errorEstimator->number_h_refinements == physics._numberHRefinements
    );
  BOOST_REQUIRE(
    physics._errorEstimator->number_p_refinements == physics._numberPRefinements
    );

}

BOOST_AUTO_TEST_CASE( PhysicsLibmesh_utils )
{
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];
  char* name = new char[27];
  strcpy(name, "dummy");
  av[0] = name;
  MPI_Init(&ac,&av);
  Communicator comm(MPI_COMM_WORLD);
  LibMeshInit libmesh_init(ac, av, comm.get()) ;
  delete av, name;

  {
  GetPot inputfile;
  inputfile = GetPot( );
  AGNOS::PhysicsLibmesh<T_S,T_P> physics(
      comm,inputfile
      );

  // we have to make a real mesh so we have some dofs
  physics._mesh = new libMesh::Mesh(physics._communicator);
  libMesh::MeshTools::Generation::build_line(
      *static_cast<libMesh::Mesh*>(physics._mesh),1,0,1,EDGE2);

  // and we need a real system to reference
  physics._equationSystems 
    = new libMesh::EquationSystems(*physics._mesh);
  physics._system = 
    &( physics._equationSystems->add_system("Basic","test") );
  physics._system->add_variable ("u", FIRST);
  physics._equationSystems->init();

  // dummy test vector to use for setting primal and adjoint solutions
  T_P testVec( std::vector<double>( 2,10.) );

  // _setPrimalSolution
  NumericVector<Number>& sol = *(physics._system->solution);
  physics._setPrimalSolution( testVec );
  BOOST_REQUIRE_CLOSE(
      sol(0), 10., 1) ;
  // _setAdjointSolution
  physics._system->add_adjoint_solution() ;
  physics._setAdjointSolution( testVec, 0 );
  
  // Testing set solVectors routine 
  // - this should also test that _set{Primal,Adjoint}Solution worked correctly
  std::map<std::string, T_P > solutionVectors ;

  physics._insertSolVec( 
      *(physics._system->solution), "primal", solutionVectors );
  BOOST_REQUIRE_CLOSE( 
      solutionVectors["primal"](0), 10., 1e-16 );

  physics._insertSolVec( 
      physics._system->get_adjoint_solution(0), "adjoint", solutionVectors );
  BOOST_REQUIRE_CLOSE( 
      solutionVectors["adjoint"](0), 10., 1e-16 );

    if( physics._mesh != NULL            ){ delete physics._mesh; }
    if( physics._meshRefinement != NULL  ){ delete physics._meshRefinement; }
    if( physics._errorEstimator != NULL  ){ delete physics._errorEstimator; }
    if( physics._qois != NULL            ){ delete physics._qois; }
    if( physics._equationSystems != NULL ){ delete physics._equationSystems; }
  }

}
