
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
#include "agnosDefines.h"
#include "Parameter.h"
#include "PseudoSpectralTensorProduct.h"
#include "PhysicsChannelFlow.h"

// boost
#define BOOST_TEST_MODULE channelTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace AGNOS;

BOOST_AUTO_TEST_CASE( channel_constructor )
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
  PhysicsChannelFlow<T_S,T_P> flowSolver(
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
  MPI_Init(&ac,&av);
  Communicator comm(MPI_COMM_WORLD);
  PETSC_COMM_WORLD = MPI_COMM_WORLD ;
  int ierr = PetscInitialize(&ac, const_cast<char***>(&av),NULL,NULL);
  LibMeshInit init (ac, av);

  GetPot inputfile;
  inputfile = GetPot( );
  inputfile.set("channel_input","flow.in") ;
  PhysicsChannelFlow<T_S,T_P> flowSolver(
      comm,inputfile
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

  BOOST_REQUIRE_CLOSE( turb->_kap, params[0]*0.41 ,1e-9);
  BOOST_REQUIRE_CLOSE( turb->_cb1 , params[1]*0.1355, 1e-9);
  BOOST_REQUIRE_CLOSE( turb->_sig  , params[2]*2.0/3.0, 1e-9);
  BOOST_REQUIRE_CLOSE( turb->_cb2    , params[3]*0.622, 1e-9);
  BOOST_REQUIRE_CLOSE( turb->_cv1 , params[4]*7.1, 1e-9);
  BOOST_REQUIRE_CLOSE( turb->_cw2 , params[5]*0.3, 1e-9);
  BOOST_REQUIRE_CLOSE( turb->_cw3 , params[6]*2.0, 1e-9);

  // make sure nothing catastrophic happens if we try to solve
  flowSolver._solve() ;

}

BOOST_AUTO_TEST_CASE( channel_convergence )
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
  inputfile.set("channel_input","flow.in") ;
  std::shared_ptr<PhysicsLibmesh<T_S,T_P> > flowSolver ; 
  flowSolver = std::shared_ptr<PhysicsChannelFlow<T_S,T_P> >(
      new PhysicsChannelFlow<T_S,T_P> ( comm,inputfile )
      );

  // construct parameter object
  unsigned int dimension = 7;
  std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;

  parameters.reserve(dimension);
  for(unsigned int i =0; i<dimension; i++)
    parameters.push_back( 
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter("CONSTANT",1.0, 1.0) )
      ); 

  // build the constantSurrogate model
  std::vector<unsigned int> order(dimension,0);
  order[0] = 4;
  order[4] = 4;
  std::set<std::string> computeSolutions ;
  computeSolutions.insert("primal");
  computeSolutions.insert("adjoint");
  /* computeSolutions.insert("errorEstimate"); */
  PseudoSpectralTensorProduct<T_S,T_P> constantSurrogate(
        comm,
        flowSolver,
        parameters, 
        order,
        computeSolutions
        );

  constantSurrogate.build( );

  // get the nominal solution we will compare to
  std::vector<double> nominalParamValue(dimension,1.);
  T_S paramValue(nominalParamValue);
  T_P primalSol = constantSurrogate.evaluate( "primal", paramValue );
  T_P adjointSol = constantSurrogate.evaluate( "adjoint", paramValue );


  // update parameters to small range around nominal values 
  parameters.clear();
  parameters.reserve(dimension);
  parameters.push_back( 
      std::shared_ptr<AGNOS::Parameter>(
        new AGNOS::Parameter("UNIFORM",0.90, 1.05) ) 
      );
  for(unsigned int i =1; i<dimension; i++)
    parameters.push_back( 
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter("CONSTANT",1.00, 1.00) )
      ); 
  parameters[4] = std::shared_ptr<AGNOS::Parameter>( new AGNOS::Parameter("UNIFORM",0.90, 1.05) )  ;
  
  // build the uniformSurrogate model
  /* std::vector<unsigned int> higherOrder(dimension,1); */
  std::shared_ptr<AGNOS::PseudoSpectralTensorProduct<T_S,T_P> >
    uniformSurrogate( new PseudoSpectralTensorProduct<T_S,T_P>(
        comm,
        flowSolver,
        parameters, 
        order  )
  );

  // build a secondary surrogate just to test errorEstimate
  std::vector<unsigned int> increaseOrder(dimension,0);
  unsigned int multiplyOrder(2);
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
  int iter=0, maxIter=4;
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

    surrogateError = secondarySurrogate->l2Norm("errorEstimate")(0);
    std::cout << "surrogateError = " << surrogateError << std::endl;

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
