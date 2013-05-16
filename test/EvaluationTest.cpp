

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Evaluation
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

  // mixed-order 5 dimensional exammple
  T_P myFunction (const T_S& paramVec)
  {
    T_P returnVec(1);
    returnVec(0) = paramVec(0) * paramVec(4) + 
      paramVec(3)*paramVec(2)*paramVec(1) + std::pow(paramVec(2),2.0);
    return returnVec ;
  }

BOOST_AUTO_TEST_SUITE(Evaluation_tensorProduct)


BOOST_AUTO_TEST_CASE(mixedPoly)
{
  BOOST_TEST_MESSAGE(" Testing evaluation of mixed order function");

  unsigned int dimension = 5;
  std::vector<unsigned int> myOrder(dimension,0);
  myOrder[0] = 1;
  myOrder[1] = 1;
  myOrder[2] = 2;
  myOrder[3] = 1;
  myOrder[4] = 1;

  std::vector<Parameter*> myParameters(
      dimension, 
      new Parameter(UNIFORM, 0.0,1.0)
      ); 

  PhysicsFunction<T_S,T_P>* myPhysicsFunction =
    new PhysicsFunctionSimple<T_S,T_P>( "myFunction", &myFunction ) ;

  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
    PseudoSpectralTensorProduct<T_S,T_P>(
        myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );

  T_S testValue(dimension);
  testValue(0) = 0.15;
  testValue(1) = 0.25;
  testValue(2) = 0.35;
  testValue(3) = 0.45;
  testValue(4) = 0.55;

  T_P surrogateValue = mySurrogate->evaluate( "myFunction", testValue);

  T_P trueValue = myFunction(testValue);

  BOOST_CHECK_CLOSE( surrogateValue(0), trueValue(0) , 1e-9 );
}

BOOST_AUTO_TEST_SUITE_END()

//________________________________________________________________//

// EOF
