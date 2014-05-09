
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

// local includes
#include <iostream>
#include <stdio.h>
#include <assert.h>

#include "agnosDefines.h"
/* #undef AGNOS_DEBUG */
/* #define AGNOS_DEBUG 0 */
#include "Parameter.h"
#include "PseudoSpectralTensorProduct.h"

#include "PhysicsCatenaryLibmesh.h"

using namespace AGNOS;

//________________________________________________________________//

  typedef libMesh::DenseVector<double> T_P ;
  typedef libMesh::DenseVector<double> T_S ;
  // linear test function
  //


  class CatenaryTest : public CppUnit::TestFixture{

    CPPUNIT_TEST_SUITE( CatenaryTest );
    CPPUNIT_TEST( runTest );
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

      void runTest( )
      {

        std::cout << "comm.size(): " << comm.size() << std::endl;
        inputfile.set("solutions","primal adjoint qoi exactQoi") ;


        unsigned int dimension = 1;
        
        std::vector<std::shared_ptr<AGNOS::Parameter> > myParameters(
            dimension, 
            std::shared_ptr<AGNOS::Parameter>(
              new AGNOS::Parameter(UNIFORM, 1.0,3.0)
            )
            ); 


        std::shared_ptr<PhysicsCatenaryLibmesh<T_S,T_P> > myPhysics( 
            new PhysicsCatenaryLibmesh<T_S,T_P>( physicsComm, inputfile ) 
            );

        std::vector<unsigned int> myOrder(dimension,4);
        PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
          new PseudoSpectralTensorProduct<T_S,T_P>(
              comm,
              myPhysics, 
              myParameters, 
              myOrder  
              );


        T_S paramValue(dimension);
        paramValue(0) = 1.5;
        T_P testValue, exactValue; 

        mySurrogate->build();
        unsigned int maxIter = 9;
        for (unsigned int iter=0; iter < maxIter; iter++)
        {
          testValue = mySurrogate->evaluate( "qoi", paramValue ) ;
          exactValue = mySurrogate->evaluate( "exactQoi", paramValue ) ;
          std::cout << "testValue.size():" << testValue.size() << std::endl;
          std::cout << "testValue = " << testValue(0) << std::endl;
          std::cout << " n_elem = " << myPhysics->getMesh( ).n_active_elem() << std::endl;

          /* myPhysics->refine( ); */
          mySurrogate->refineUniformly();
          mySurrogate->build();
        }
        std::cout << "testValue = " << testValue(0) << std::endl;
        std::cout << " n_elem = " << myPhysics->getMesh( ).n_active_elem() << std::endl;


        CPPUNIT_ASSERT( 
            std::abs(testValue(0) - exactValue(0) ) <= 1e-6 
            );

        delete mySurrogate;
      }

  };

    CPPUNIT_TEST_SUITE_REGISTRATION( CatenaryTest );


