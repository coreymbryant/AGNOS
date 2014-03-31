
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
  std::set<std::string> computeSolutions ;
  computeSolutions.insert("primal");
  computeSolutions.insert("adjoint");
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
  for(unsigned int i =0; i<dimension; i++)
    parameters.push_back( 
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter("UNIFORM",0.90, 1.05) )
      ); 
  
  // build the uniformSurrogate model
  /* std::vector<unsigned int> higherOrder(dimension,1); */
  PseudoSpectralTensorProduct<T_S,T_P> uniformSurrogate(
        comm,
        flowSolver,
        parameters, 
        order  
        );

  T_P primalPred, adjointPred;

  int iter=0, maxIter=2;
  double tol = 1e-3, primalDiff,adjointDiff;
  do 
  {
    iter++;
    uniformSurrogate.refine();
    uniformSurrogate.build( );
    primalPred = uniformSurrogate.evaluate( "primal", paramValue );
    primalPred -= primalSol ;
    primalDiff = primalPred.linfty_norm() ;
    std::cout << "primalDiff = " << primalDiff << std::endl;
    adjointPred = uniformSurrogate.evaluate( "adjoint", paramValue );
    adjointPred -= adjointSol ;
    adjointDiff = adjointPred.linfty_norm() ;
    std::cout << "adjointDiff = " << adjointDiff << std::endl;
  } while ( 
       (iter<maxIter) && 
       (primalDiff > tol) &&
       (adjointDiff > tol) 
       ) ;

  // for 2 iterations diff should be within 1e-3
  // for 3 its within 1e-6 (this takes a very long time to run though)
  BOOST_REQUIRE_SMALL( primalDiff, tol  );
  BOOST_REQUIRE_SMALL( adjointDiff, tol  );

}
