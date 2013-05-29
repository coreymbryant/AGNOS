

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Catenary
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// local includes
#include <iostream>
#include <stdio.h>
#include <assert.h>

#include "agnosDefines.h"
#include "Parameter.h"
#include "PseudoSpectralTensorProduct.h"
#include "Parameter.h"

#include "PhysicsCatenary.h"

using namespace AGNOS;

//________________________________________________________________//

  typedef libMesh::DenseVector<double> T_P ;
  typedef libMesh::DenseVector<double> T_S ;
  // linear test function

BOOST_AUTO_TEST_SUITE(Catenary_tensorProduct)

  const Communicator comm( MPI_COMM_NULL );

  unsigned int dimension = 1;
  
  std::vector<Parameter*> myParameters(
      dimension, 
      new Parameter(UNIFORM, 1.0,3.0)
      ); 

  PhysicsModel<T_S,T_P>* myPhysics = new PhysicsCatenary<T_S,T_P>( -10.0 );
  PhysicsFunction<T_S,T_P>* myPhysicsFunction =
    new PhysicsFunctionPrimal<T_S,T_P>( *myPhysics ) ;

BOOST_AUTO_TEST_CASE(Catenary_N0)
{
  BOOST_TEST_MESSAGE(" Testing Catenary with N=0");

  std::vector<unsigned int> myOrder(dimension,0);
  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
    new PseudoSpectralTensorProduct<T_S,T_P>(
        &comm,
        myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::map< std::string, std::vector<T_P> > myCoeff 
    = mySurrogate->getCoefficients( );

  BOOST_CHECK_CLOSE( myCoeff["primal"][0](0) , -10.0/16.0, 1e-9 );

}

BOOST_AUTO_TEST_CASE(Catenary_N1)
{
  BOOST_TEST_MESSAGE(" Testing Catenary with N=1");

  std::vector<unsigned int> myOrder(dimension,1);
  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
    new PseudoSpectralTensorProduct<T_S,T_P>(
        &comm,
        myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::map< std::string, std::vector<T_P> > myCoeff 
    = mySurrogate->getCoefficients( );


  // 12/11 is the sum of u(\xi_j) 
  BOOST_CHECK_CLOSE( myCoeff["primal"][0](0) , -10.0/16.0 * (12./11.), 1e-9 );
  // -2 sqrt(3)/11 is correct contribution from poly evals
  BOOST_CHECK_CLOSE( 
      myCoeff["primal"][1](0) , 
      -10.0/16.0 * ( -2.0 * std::sqrt(3.0) / 11.0 ), 
      1e-9 );

}

BOOST_AUTO_TEST_CASE(Catenary_N4)
{
  BOOST_TEST_MESSAGE(" Testing Catenary with N=4");

  std::vector<unsigned int> myOrder(dimension,4);
  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
    new PseudoSpectralTensorProduct<T_S,T_P>(
        &comm,
        myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::map< std::string, std::vector<T_P> > myCoeff 
    = mySurrogate->getCoefficients( );


  // Coefficients generated from pmpack for comparison
  BOOST_CHECK_CLOSE( 
      myCoeff["primal"][0](0) , 
      -6.866307761327950e-01,
      1e-9 );
  BOOST_CHECK_CLOSE( 
      myCoeff["primal"][1](0) , 
      2.134952711438086e-01,
      1e-9 );
  BOOST_CHECK_CLOSE( 
      myCoeff["primal"][2](0) , 
      -5.918708419582225e-02, 
      1e-9 );
  BOOST_CHECK_CLOSE( 
      myCoeff["primal"][3](0) , 
      1.602406581398491e-02,
      1e-9 );
  BOOST_CHECK_CLOSE( 
      myCoeff["primal"][4](0) , 
      -4.037685060565489e-03,
      1e-9 );

}

BOOST_AUTO_TEST_CASE(Catenary_convergence)
{
  BOOST_TEST_MESSAGE(" Testing Catenary convergence");

  std::vector<unsigned int> myOrder(dimension,0);
  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
    new PseudoSpectralTensorProduct<T_S,T_P>(
        &comm,
        myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  T_S paramValue(dimension);
  paramValue(0) = 1.5;
  T_P testValue; 

  unsigned int maxIter = 25;
  for (unsigned int iter=0; iter < maxIter-1; iter++)
  {
    mySurrogate->refine();
    testValue = mySurrogate->evaluate( "primal", paramValue ) ;
  }

  BOOST_CHECK_CLOSE(  testValue(0) , -10.0/(8.0 * paramValue(0) ), 1e-9 );

}

BOOST_AUTO_TEST_SUITE_END()

//________________________________________________________________//

// EOF
