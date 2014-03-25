
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

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

  class EvaluationTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( EvaluationTest );
    CPPUNIT_TEST( mixedPoly );
    CPPUNIT_TEST_SUITE_END();

    public:
      Communicator comm, physicsComm;
      GetPot inputfile;

      void setUp( )
      {
        comm = Communicator(MPI_COMM_WORLD);
        inputfile = GetPot() ;

        MPI_Comm subComm;
        int mpiSplit =  
          MPI_Comm_split( MPI_COMM_WORLD, comm.rank(), 0, &subComm);
        physicsComm = Communicator(subComm) ;
      }

      void tearDown( )
      {
      }


      void mixedPoly()
      {
        std::cout << " Testing evaluation of mixed order function"<< std::endl;

        unsigned int dimension = 5;
        std::vector<unsigned int> myOrder(dimension,0);
        myOrder[0] = 1;
        myOrder[1] = 1;
        myOrder[2] = 2;
        myOrder[3] = 1;
        myOrder[4] = 1;

        std::vector<std::shared_ptr<AGNOS::Parameter> > myParameters(
            dimension, 
            std::shared_ptr<AGNOS::Parameter>(
              new AGNOS::Parameter(UNIFORM, 0.0,1.0)
              )
            ); 

        std::shared_ptr<PhysicsUser<T_S,T_P> > myPhysics(
          new PhysicsUser<T_S,T_P>(physicsComm,inputfile));
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

        CPPUNIT_ASSERT( std::abs( 
              surrogateValue(0) - trueSolution["primal"](0) ) <= 1e-9 );

      delete mySurrogate;
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( EvaluationTest );

//________________________________________________________________//

// EOF
