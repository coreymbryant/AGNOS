

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Polynomials
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
  T_P linearFunction (const T_S& paramVec)
  {
    T_P returnVec(1);
    returnVec(0)=paramVec(0);
    for(unsigned int dim=1; dim<paramVec.size(); dim++)
      returnVec(0) *= paramVec(dim) ;
    return returnVec ;
  }

  // quadratic test function
  T_P quadFunction (const T_S& paramVec)
  {
    T_P returnVec(1);
    returnVec(0) = paramVec(0) * paramVec(0) ;
    for(unsigned int dim=1; dim<paramVec.size(); dim++)
      returnVec(0) *= (paramVec(dim)*paramVec(dim)) ;
    return returnVec ;
  }

  // mixed-order 5 dimensional exammple
  T_P mixedFunction (const T_S& paramVec)
  {
    T_P returnVec(1);
    returnVec(0) = paramVec(2)  + paramVec(1)*paramVec(1)*paramVec(0) + std::pow(paramVec(3),3) ;
    return returnVec ;
  }

BOOST_AUTO_TEST_SUITE(Polynomials_tensorProduct)


BOOST_AUTO_TEST_CASE(Linear_1D)
{
  BOOST_TEST_MESSAGE(" Testing 1D linear approximation");

  unsigned int dimension = 1;
  std::vector<unsigned int> myOrder(dimension,3);

  std::vector<Parameter*> myParameters(
      dimension, 
      new Parameter(UNIFORM, -1.0,1.0)
      ); 

  PhysicsFunction<T_S,T_P>* myPhysicsFunction =
    new PhysicsFunctionSimple<T_S,T_P>( &linearFunction ) ;

  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
    PseudoSpectralTensorProduct<T_S,T_P>(
        *myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::vector<T_P> myCoeff = mySurrogate->getCoefficients( );

  BOOST_CHECK_CLOSE( myCoeff[1](0) , std::sqrt(1.0/3.0) , 1e-9 );
}

BOOST_AUTO_TEST_CASE(Linear_ND)
{
  BOOST_TEST_MESSAGE(" Testing 7D linear approximation");

  unsigned int dimension = 7;
  std::vector<unsigned int> myOrder(dimension,1);
  std::vector<Parameter*> myParameters(
      dimension, 
      new Parameter(UNIFORM, -1.0,1.0)
      ); 

  PhysicsFunction<T_S,T_P>* myPhysicsFunction =
    new PhysicsFunctionSimple<T_S,T_P>( &linearFunction ) ;

  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
    PseudoSpectralTensorProduct<T_S,T_P>(
        *myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::vector<T_P> myCoeff = mySurrogate->getCoefficients( );


  BOOST_CHECK_CLOSE( 
      myCoeff.back()(0) , 
      std::pow(std::sqrt(1.0/3.0),dimension),
      1e-4 );
}

BOOST_AUTO_TEST_CASE(Quad_1D)
{
  BOOST_TEST_MESSAGE(" Testing 1D quadratic approximation");

  unsigned int dimension = 1;
  std::vector<unsigned int> myOrder(dimension,2);

  std::vector<Parameter*> myParameters(
      dimension, 
      new Parameter(UNIFORM, -1.0,1.0)
      ); 

  PhysicsFunction<T_S,T_P>* myPhysicsFunction =
    new PhysicsFunctionSimple<T_S,T_P>( &quadFunction ) ;

  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
    PseudoSpectralTensorProduct<T_S,T_P>(
        *myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::vector<T_P> myCoeff = mySurrogate->getCoefficients( );

  BOOST_CHECK_CLOSE( myCoeff[0](0) , 1.0/3.0 , 1e-9 );
  BOOST_CHECK_CLOSE( myCoeff[2](0) , 2.0/3.0 * std::sqrt(1.0/5.0), 1e-9 );
}

BOOST_AUTO_TEST_CASE(Quad_ND)
{
  BOOST_TEST_MESSAGE(" Testing 3D quadratic approximation");

  unsigned int dimension = 3;
  std::vector<unsigned int> myOrder(dimension,2);

  std::vector<Parameter*> myParameters(
      dimension, 
      new Parameter(UNIFORM, -1.0,1.0)
      ); 

  PhysicsFunction<T_S,T_P>* myPhysicsFunction =
    new PhysicsFunctionSimple<T_S,T_P>( &quadFunction ) ;

  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
    PseudoSpectralTensorProduct<T_S,T_P>(
        *myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::vector<T_P> myCoeff = mySurrogate->getCoefficients( );

  BOOST_CHECK_CLOSE( myCoeff[0](0) , std::pow(1.0/3.0,dimension), 1e-9 );
  BOOST_CHECK_CLOSE( 
      myCoeff.back()(0) , 
      std::pow( std::sqrt(1.0/5.0) * 2.0/3.0 ,dimension), 
      1e-9 );
}

BOOST_AUTO_TEST_CASE(OrderN_1D)
{
  BOOST_TEST_MESSAGE(" Testing mixed-order 5 dimensional example");

  unsigned int dimension = 5;
  std::vector<unsigned int> myOrder(dimension,0);
  myOrder[0] = 1;
  myOrder[1] = 2;
  myOrder[2] = 1;
  myOrder[3] = 3;
  myOrder[4] = 0;

  std::vector<Parameter*> myParameters(
      dimension, 
      new Parameter(UNIFORM, -1.0,1.0)
      ); 

  PhysicsFunction<T_S,T_P>* myPhysicsFunction =
    new PhysicsFunctionSimple<T_S,T_P>( &mixedFunction ) ;

  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
    PseudoSpectralTensorProduct<T_S,T_P>(
        *myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::vector<T_P> myCoeff = mySurrogate->getCoefficients( );

  /* mySurrogate->printIntegrationWeights() ; */
  /* mySurrogate->printIntegrationPoints() ; */
  /* mySurrogate->printIndexSet() ; */
  /* for(unsigned int coeff=0; coeff<myCoeff.size(); coeff++) */
  /*   for(unsigned int dim=0; dim<myCoeff[coeff].size(); dim++) */
  /*     std::cout << "coeff[" << coeff << "](" << dim << ") = " */ 
  /*       << myCoeff[coeff](dim) << std::endl; */

  // 0 0 1 0 0  = 1 * 1/norm(psi) from normailziation
  BOOST_CHECK_CLOSE( 
      myCoeff[4](0) , 
      std::sqrt(1.0/3.0),
      1e-9 );

  // 1 2 0 0 0  = 1
  // have to use combination of terms since psi_2 = 1/2*(3x^2-1) 
  BOOST_CHECK_CLOSE(
      myCoeff[40](0),
      2.0/(3.0 * std::sqrt(3.0*5.0)) ,
      1e-9 );
  BOOST_CHECK_CLOSE(
      myCoeff[24](0),
      1.0 / ( 3.0 * std::sqrt(3.0) ),
      1e-9 );

  // 0 0 0 3 0  = 1
  // have to use combination of terms since psi_2 = 1/2*(5x^3-3x) 
  BOOST_CHECK_CLOSE(
      myCoeff[3](0),
      2.0 / ( 5.0 * std::sqrt(7.0) ),
      1e-9 );
  BOOST_CHECK_CLOSE(
       myCoeff[1](0),
       3.0 / ( 5.0 * std::sqrt(3.0) ), 
       1e-9 );
}

BOOST_AUTO_TEST_SUITE_END()

//________________________________________________________________//

// EOF
