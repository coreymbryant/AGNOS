

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Catenary
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// local includes
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include "Parameter.h"
#include "PseudoSpectralTensorProduct.h"
#include "PhysicsCatenary.h"
#include "PhysicsFunction.h"

#include "libmesh/dense_vector.h"
using namespace AGNOS;

//________________________________________________________________//

  typedef libMesh::DenseVector<double> T_P ;
  typedef libMesh::DenseVector<double> T_S ;
  // linear test function

BOOST_AUTO_TEST_SUITE(Catenary_tensorProduct)


  unsigned int dimension = 1;
  
  std::vector<Parameter*> myParameters(
      dimension, 
      new Parameter(UNIFORM, 1.0,3.0)
      ); 

  PhysicsModel<T_S,T_P>* myPhysics = new PhysicsCatenary<T_S,T_P>( );
  PhysicsFunction<T_S,T_P>* myPhysicsFunction =
    new PhysicsFunctionPrimal<T_S,T_P>( *myPhysics ) ;

BOOST_AUTO_TEST_CASE(Catenary_N0)
{
  BOOST_TEST_MESSAGE(" Testing Catenary with N=0");

  std::vector<unsigned int> myOrder(dimension,0);
  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
    new PseudoSpectralTensorProduct<T_S,T_P>(
        *myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::vector<T_P> myCoeff = mySurrogate->getCoefficients( );

  BOOST_CHECK_CLOSE( myCoeff[0](0) , -10.0/16.0, 1e-9 );

}

BOOST_AUTO_TEST_CASE(Catenary_N1)
{
  BOOST_TEST_MESSAGE(" Testing Catenary with N=1");

  std::vector<unsigned int> myOrder(dimension,1);
  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
    new PseudoSpectralTensorProduct<T_S,T_P>(
        *myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::vector<T_P> myCoeff = mySurrogate->getCoefficients( );


  // 12/11 is the sum of u(\xi_j) 
  BOOST_CHECK_CLOSE( myCoeff[0](0) , -10.0/16.0 * (12./11.), 1e-9 );
  // -2 sqrt(3)/11 is correct contribution from poly evals
  BOOST_CHECK_CLOSE( 
      myCoeff[1](0) , 
      -10.0/16.0 * ( -2.0 * std::sqrt(3.0) / 11.0 ), 
      1e-9 );

}

BOOST_AUTO_TEST_CASE(Catenary_N5)
{
  BOOST_TEST_MESSAGE(" Testing Catenary with N=4");

  std::vector<unsigned int> myOrder(dimension,4);
  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
    new PseudoSpectralTensorProduct<T_S,T_P>(
        *myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::vector<T_P> myCoeff = mySurrogate->getCoefficients( );

  /* mySurrogate->printIntegrationPoints(); */
  /* mySurrogate->printIntegrationWeights(); */
  /* for(unsigned int coeff=0; coeff<myCoeff.size(); coeff++) */
  /*     std::cout << std::setprecision(5) << std::scientific */ 
  /*       << "coeff[" << coeff << "](0) = " */ 
  /*       << myCoeff[coeff](0) << std::endl; */

  // Coefficients generated from pmpack for comparison
  BOOST_CHECK_CLOSE( 
      myCoeff[0](0) , 
      -6.866307761327950e-01,
      1e-9 );
  BOOST_CHECK_CLOSE( 
      myCoeff[1](0) , 
      2.134952711438086e-01,
      1e-9 );
  BOOST_CHECK_CLOSE( 
      myCoeff[2](0) , 
      -5.918708419582225e-02, 
      1e-9 );
  BOOST_CHECK_CLOSE( 
      myCoeff[3](0) , 
      1.602406581398491e-02,
      1e-9 );
  BOOST_CHECK_CLOSE( 
      myCoeff[4](0) , 
      -4.037685060565489e-03,
      1e-9 );

}

BOOST_AUTO_TEST_SUITE_END()

//________________________________________________________________//

// EOF
