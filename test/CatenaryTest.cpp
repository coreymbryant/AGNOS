

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

// local includes
#include <iostream>
#include <stdio.h>
#include <assert.h>

#include "agnosDefines.h"
#undef AGNOS_DEBUG
#define AGNOS_DEBUG 0


#include "Parameter.h"
#include "PseudoSpectralTensorProduct.h"

#include "PhysicsCatenary.h"


using namespace AGNOS;

//________________________________________________________________//

  typedef libMesh::DenseVector<double> T_P ;
  typedef libMesh::DenseVector<double> T_S ;
  // linear test function

  class CatenaryTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( CatenaryTest );
    CPPUNIT_TEST( catenaryN0 );
    CPPUNIT_TEST( catenaryN1 );
    CPPUNIT_TEST( catenaryN4 );
    CPPUNIT_TEST( catenaryConvergence );
    CPPUNIT_TEST_SUITE_END();

    public:
      Communicator comm, physicsComm;
      GetPot inputfile;
      std::shared_ptr<PhysicsModel<T_S,T_P> > myPhysics ; 
      unsigned int dimension ;
        
      std::vector<std::shared_ptr<AGNOS::Parameter> > myParameters;

      void setUp( )
      {
        comm = Communicator(MPI_COMM_WORLD);
        inputfile = GetPot() ;

        MPI_Comm subComm;
        int mpiSplit =  
          MPI_Comm_split( MPI_COMM_WORLD, comm.rank(), 0, &subComm);
        physicsComm = Communicator(subComm) ;

        dimension = 1;
        myPhysics = std::shared_ptr<PhysicsModel<T_S,T_P> > (
            new PhysicsCatenary<T_S,T_P>(physicsComm,inputfile ) 
            );
        myParameters.reserve(1);
        myParameters.push_back( 
            std::shared_ptr<AGNOS::Parameter>(
              new AGNOS::Parameter(UNIFORM, 1.0,3.0)
            )
          ); 
      }

      void tearDown( )
      {
      }

      
      void catenaryN0()
      {
        std::cout << "\n--------------------------\n Testing Catenary with N=0"
          << std::endl;

        std::vector<unsigned int> myOrder(dimension,0);
        PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
          new PseudoSpectralTensorProduct<T_S,T_P>(
              comm,
              myPhysics,
              myParameters, 
              myOrder  
              );

        mySurrogate->build( );
        std::map< std::string, LocalMatrix > myCoeff 
          = mySurrogate->getCoefficients( );


        CPPUNIT_ASSERT( 
            std::abs(myCoeff["primal"](0,0)  - -10.0/16.0) <= 1e-9
            );

        delete mySurrogate;
      }



      void catenaryN1()
      {
        std::cout << "\n--------------------------\n Testing Catenary with N=1"
          << std::endl;

        std::vector<unsigned int> myOrder(dimension,1);
        PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
          new PseudoSpectralTensorProduct<T_S,T_P>(
              comm,
              myPhysics, 
              myParameters, 
              myOrder  
              );

        mySurrogate->build( );
        std::map< std::string, LocalMatrix > myCoeff 
          = mySurrogate->getCoefficients( );


        // 12/11 is the sum of u(\xi_j) 
        CPPUNIT_ASSERT( 
            std::abs(myCoeff["primal"](0,0)  - -10.0/16.0 * (12./11.)) <= 1e-9 
            );
        // -2 sqrt(3)/11 is correct contribution from poly evals
        CPPUNIT_ASSERT( 
            std::abs(
              myCoeff["primal"](1,0)  - 
              -10.0/16.0 * ( -2.0 * std::sqrt(3.0) / 11.0 ) 
              ) <= 1e-9 );

        delete mySurrogate;
      }

      void catenaryN4()
      {
        std::cout << "\n--------------------------\n Testing Catenary with N=4"
          << std::endl;

        std::vector<unsigned int> myOrder(dimension,4);
        PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
          new PseudoSpectralTensorProduct<T_S,T_P>(
              comm,
              myPhysics, 
              myParameters, 
              myOrder  
              );

        mySurrogate->build( );
        std::map< std::string, LocalMatrix > myCoeff 
          = mySurrogate->getCoefficients( );


        // Coefficients generated from pmpack for comparison
        CPPUNIT_ASSERT( 
            std::abs(myCoeff["primal"](0,0)  - 
            -6.866307761327950e-01 )
            <= 1e-9 );
        CPPUNIT_ASSERT( 
            std::abs(myCoeff["primal"](1,0)  - 
            2.134952711438086e-01 ) <=
            1e-9 );
        CPPUNIT_ASSERT( 
            std::abs(myCoeff["primal"](2,0)  - 
            -5.918708419582225e-02 ) <= 
            1e-9 );
        CPPUNIT_ASSERT( 
            std::abs( myCoeff["primal"](3,0)  - 
            1.602406581398491e-02 ) <=
            1e-9 );
        CPPUNIT_ASSERT( 
            std::abs( myCoeff["primal"](4,0)  - 
            -4.037685060565489e-03 ) <=
            1e-9 );

        delete mySurrogate;
      }

      void catenaryConvergence()
      {
        std::cout 
          << "\n--------------------------\n Testing Catenary convergence" 
          << std::endl;

        std::vector<unsigned int> myOrder(dimension,0);
        PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
          new PseudoSpectralTensorProduct<T_S,T_P>(
              comm,
              myPhysics, 
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
          mySurrogate->build();
          testValue = mySurrogate->evaluate( "primal", paramValue ) ;
        }

        /* std::cout << "eval = " << testValue(0) << std::endl; */
        CPPUNIT_ASSERT(  std::abs(testValue(0) - -10.0/(8.0 * paramValue(0) )) 
            <=  1e-9 );

        delete mySurrogate;
      }


  };

  CPPUNIT_TEST_SUITE_REGISTRATION( CatenaryTest );

//________________________________________________________________//

// EOF
