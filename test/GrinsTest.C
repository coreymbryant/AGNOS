
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
#include "agnosDefines.h"
#include "Parameter.h"
#include "PseudoSpectralTensorProduct.h"
#include "PhysicsGrins.h"
#include "grins/inc_navier_stokes.h"

// boost
#define BOOST_TEST_MODULE grinsTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace AGNOS;

BOOST_AUTO_TEST_CASE( grins_constructor )
{
  
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];

  char* name = new char[27];
  strcpy(name, "dummy");
  
  av[0] = name;
  
  // Initialize libmesh
  MPI_Init(&ac,&av);
  Communicator comm(MPI_COMM_WORLD);
  PETSC_COMM_WORLD = MPI_COMM_WORLD ;
  int ierr = PetscInitialize(&ac, const_cast<char***>(&av),NULL,NULL);
  LibMeshInit init (ac, av);

  GetPot inputfile;
  inputfile = GetPot( );
  PhysicsGrins<T_S,T_P> grinsSolver( comm, inputfile );

  BOOST_REQUIRE( grinsSolver._availableSolutions.count("primal") );
  BOOST_REQUIRE( grinsSolver._availableSolutions.count("adjoint") );
  BOOST_REQUIRE( grinsSolver._availableSolutions.count("qoi") );
  BOOST_REQUIRE( grinsSolver._availableSolutions.count("errorEstimate") );
  BOOST_REQUIRE( grinsSolver._availableSolutions.count("errorIndicators") );
  
  // check for initial parameters set to NULL values
  BOOST_REQUIRE( &grinsSolver.getSystem() != NULL );
  BOOST_REQUIRE( &grinsSolver.getEquationSystems() != NULL );
  BOOST_REQUIRE( &grinsSolver.getMesh() != NULL );
  BOOST_REQUIRE( &grinsSolver.getMeshRefinement() != NULL );
  BOOST_REQUIRE( &grinsSolver.getEstimator() != NULL );
  BOOST_REQUIRE( &grinsSolver.getQois() != NULL );
}

BOOST_AUTO_TEST_CASE( grins_solve )
{
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];

  char* name = new char[27];
  strcpy(name, "dummy");
  
  av[0] = name;
  
  // Initialize libmesh
  MPI_Init(&ac,&av);
  Communicator comm(MPI_COMM_WORLD);
  PETSC_COMM_WORLD = MPI_COMM_WORLD ;
  int ierr = PetscInitialize(&ac, const_cast<char***>(&av),NULL,NULL);
  LibMeshInit init (ac, av);

  GetPot inputfile;
  inputfile = GetPot( );
  inputfile.set("physics/grins_input","./grins.in") ;
  inputfile.set("physics/IncompressibleNavierStokes/rho","$(0)") ;
  inputfile.set("physics/IncompressibleNavierStokes/mu","$(1)") ;
  PhysicsGrins<T_S,T_P> grinsSolver(
      comm,inputfile
      );

  // test setting parameter values
  std::vector<double> params(2,1.);
  T_S parameterValues(params) ;
  grinsSolver._setParameterValues(parameterValues);
  GRINS::IncompressibleNavierStokes* physics = (
      static_cast<GRINS::IncompressibleNavierStokes*>(
        grinsSolver._multiphysicsSystem
        ->get_physics("IncompressibleNavierStokes").get() )
      ) ;

  std::cout <<" rho: " << physics->_rho << std::endl;
  std::cout <<" mu: " << physics->_mu << std::endl;
  BOOST_REQUIRE_CLOSE( physics->_rho, 1.0 ,1e-9);
  BOOST_REQUIRE_CLOSE( physics->_mu , 1.0, 1e-9);

  // make sure nothing catastrophic happens if we try to solve
  grinsSolver._solve() ;

}

BOOST_AUTO_TEST_CASE( grins_convergence )
{
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];

  char* name = new char[27];
  strcpy(name, "dummy");
  
  av[0] = name;
  
  // Initialize libmesh
  MPI_Init(&ac,&av);
  Communicator comm(MPI_COMM_WORLD);
  PETSC_COMM_WORLD = MPI_COMM_WORLD ;
  int ierr = PetscInitialize(&ac, const_cast<char***>(&av),NULL,NULL);
  LibMeshInit init (ac, av);

  // istantiate phyiscs model
  GetPot inputfile;
  inputfile = GetPot( );
  inputfile.set("physics/grins_input","./grins.in") ;
  /* inputfile.set("physics/IncompressibleNavierStokes/mu","$(0)") ; */
  /* inputfile.set( */
  /*     "physics/IncompressibleNavierStokes/parabolic_profile_coeffs_1", */
  /*     "0.0 0.0 -$(1)/32.0 0.0 0.0 $(1)/8.0") ; */

  std::shared_ptr<PhysicsLibmesh<T_S,T_P> > grinsSolver ; 
  grinsSolver = std::shared_ptr<PhysicsGrins<T_S,T_P> >(
      new PhysicsGrins<T_S,T_P> ( comm,inputfile )
      );

  // construct parameter object
  unsigned int dimension = 2;
  std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;


  // build the constantSurrogate model
  std::vector<unsigned int> order(dimension,2);
  parameters.push_back( 
      std::shared_ptr<AGNOS::Parameter>(
        new AGNOS::Parameter("CONSTANT",0.055, 0.055) )
    ); 
  parameters.push_back( 
      std::shared_ptr<AGNOS::Parameter>(
        new AGNOS::Parameter("CONSTANT",6.0, 6.0) )
    ); 

  std::set<std::string> computeSolutions ;
  computeSolutions.insert("primal");
  computeSolutions.insert("adjoint");
  computeSolutions.insert("errorEstimate");
  PseudoSpectralTensorProduct<T_S,T_P> constantSurrogate(
        comm,
        grinsSolver,
        parameters, 
        order,
        computeSolutions
        );

  constantSurrogate.build( );

  // get the nominal solution we will compare to
  std::vector<double> nominalParamValue(dimension,1.);
  // neither of these is the mid point
  nominalParamValue[0] = 0.05 ; 
  nominalParamValue[1] = 5.00 ;
  T_S paramValue(nominalParamValue);
  T_P primalSol = constantSurrogate.evaluate( "primal", paramValue );
  T_P adjointSol = constantSurrogate.evaluate( "adjoint", paramValue );


  // update parameters to small range around nominal values 
  parameters.clear();
  parameters.reserve(dimension);
  parameters.push_back( 
      std::shared_ptr<AGNOS::Parameter>(
        new AGNOS::Parameter("UNIFORM",0.01, 0.1) )
    ); 
  parameters.push_back( 
      std::shared_ptr<AGNOS::Parameter>(
        new AGNOS::Parameter("UNIFORM",3.0, 9.0) )
    ); 
  
  // build the uniformSurrogate model
  /* std::vector<unsigned int> higherOrder(dimension,1); */
  std::shared_ptr<AGNOS::PseudoSpectralTensorProduct<T_S,T_P> >
    uniformSurrogate( new PseudoSpectralTensorProduct<T_S,T_P>(
        comm,
        grinsSolver,
        parameters, 
        order  )
  );

  // build a secondary surrogate just to test errorEstimate
  std::vector<unsigned int> increaseOrder(dimension,1);
  unsigned int multiplyOrder(1);
  std::set<std::string> secondarySolutions ;
  secondarySolutions.insert("errorEstimate");
  std::shared_ptr<AGNOS::PseudoSpectralTensorProduct<T_S,T_P> >
  secondarySurrogate( new PseudoSpectralTensorProduct<T_S,T_P>(
      uniformSurrogate,
      increaseOrder,
      multiplyOrder,
      computeSolutions,
      secondarySolutions )
      );

  T_P primalPred, adjointPred;
  double surrogateError = 0.;
  int iter=0, maxIter=2;
  double tol = 1e-3, primalDiff,adjointDiff;
  bool condition = true ;
  do 
  {
    iter++;

    uniformSurrogate->refine( );
    secondarySurrogate->refine();
    uniformSurrogate->build( );
    secondarySurrogate->build( );

    primalPred = uniformSurrogate->evaluate( "primal", paramValue );
    primalPred -= primalSol ;
    primalDiff = primalPred.linfty_norm() ;
    std::cout << "primalDiff = " << primalDiff << std::endl;

    adjointPred = uniformSurrogate->evaluate( "adjoint", paramValue );
    adjointPred -= adjointSol ;
    adjointDiff = adjointPred.linfty_norm() ;
    std::cout << "adjointDiff = " << adjointDiff << std::endl;

    surrogateError = uniformSurrogate->l2NormDifference(
        *secondarySurrogate,
        "errorEstimate");
    std::cout << "physicalError = " <<
      uniformSurrogate->l2Norm("errorEstimate")(0) << std::endl;
    std::cout << "surrogateError = " << surrogateError << std::endl;
    std::cout << "totalError = " <<
      secondarySurrogate->l2Norm("errorEstimate")(0) << std::endl;

    std::cout << " iter<maxIter: " << (iter<maxIter) << std::endl;
    std::cout << " primalDiff > tol: " << (primalDiff>tol) << std::endl;
    std::cout << " adjointDiff > tol: " << (adjointDiff>tol) << std::endl;
    std::cout << " surrogateError > tol: " << (surrogateError>tol) << std::endl;
    condition = 
       (iter<maxIter) && (
       (primalDiff > tol) || (adjointDiff > tol) || (surrogateError > tol) 
       );
    std::cout << "condition: " << condition << std::endl;

  } while ( condition ) ;

  // for 2 iterations diff should be within 1e-3
  // for 3 its within 1e-6 (this takes a very long time to run though)
  BOOST_REQUIRE_SMALL( primalDiff, tol  );
  BOOST_REQUIRE_SMALL( adjointDiff, tol  );
  BOOST_REQUIRE_SMALL( surrogateError, tol  );

}
