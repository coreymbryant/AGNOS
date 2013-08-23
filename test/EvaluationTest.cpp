

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Evaluation
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

  // linear test function

  // mixed-order 5 dimensional exammple
  void myFunction (
      const T_S& paramVec,
      std::map<std::string, T_P>& solutionVectors
      )
  {
    T_P returnVec(1);
    returnVec(0) = paramVec(0) * paramVec(4) + 
      paramVec(3)*paramVec(2)*paramVec(1) + std::pow(paramVec(2),2.0);

    solutionVectors.clear();
    solutionVectors.insert(std::pair<std::string,T_P>("primal",returnVec) );
  }

BOOST_AUTO_TEST_SUITE(Evaluation_tensorProduct)

  const Communicator comm( MPI_COMM_NULL );
  const GetPot inputfile = GetPot();

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

  PhysicsModel<T_S,T_P>* myPhysics=
    new PhysicsModel<T_S,T_P>(comm,inputfile);
  myPhysics->attach_compute_function(&myFunction);

  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
    PseudoSpectralTensorProduct<T_S,T_P>(
        comm,
        myPhysics, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );

  std::map<std::string,T_P> trueSolution;
  T_S testValue(dimension);
  testValue(0) = 0.15;
  testValue(1) = 0.25;
  testValue(2) = 0.35;
  testValue(3) = 0.45;
  testValue(4) = 0.55;

  T_P surrogateValue = mySurrogate->evaluate( "primal", testValue);

  myFunction(testValue,trueSolution);

  BOOST_CHECK_CLOSE( surrogateValue(0), trueSolution["primal"](0) , 1e-9 );
}

BOOST_AUTO_TEST_SUITE_END()

//________________________________________________________________//

// EOF
