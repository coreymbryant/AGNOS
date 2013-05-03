

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Polynomials
#include <boost/test/unit_test.hpp>

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

BOOST_AUTO_TEST_SUITE(Polynomials_tensorProduct)

  typedef libMesh::DenseVector<double> T_P ;
  typedef libMesh::DenseVector<double> T_S ;
  T_P myFunction (T_S& paramVec)
  {
    T_P dummVec = paramVec ;
    return dummyVec ;
  }

BOOST_AUTO_TEST_CASE(Linear_1D)
{

  unsigned int dimension = 1;

  std::vector<Parameter*> myParameters(
      dimension, 
      new Parameter(UNIFORM, -1.0,1.0)
      ); 

  std::vector<unsigned int> myOrder(dimension,1);

  PhysicsFunction<T_S,T_P>* myPhysicsFunction =
    new PhysicsFunctionSimple<T_S,T_P>( &myFunction ) ;

  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
    PseudoSpectralTensorProduct<T_S,T_P>(
        *myPhysicsFunction, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );

  BOOST_REQUIRE( 1 );
}

BOOST_AUTO_TEST_CASE(Linear_ND)
{
  BOOST_REQUIRE( 1 );
}

BOOST_AUTO_TEST_CASE(Quad_1D)
{
  BOOST_REQUIRE( 1 );
}

BOOST_AUTO_TEST_CASE(Quad_ND)
{
  BOOST_REQUIRE( 1 );
}

BOOST_AUTO_TEST_CASE(OrderN_1D)
{
  BOOST_REQUIRE( 1 );
}

BOOST_AUTO_TEST_SUITE_END()

//________________________________________________________________//

// EOF
