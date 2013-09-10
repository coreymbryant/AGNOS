

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

#include "PhysicsViscousBurgers.h"


using namespace AGNOS;

//________________________________________________________________//

  typedef libMesh::DenseVector<double> T_P ;
  typedef libMesh::DenseVector<double> T_S ;
  // linear test function

  class ViscousBurgersTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( ViscousBurgersTest );
    CPPUNIT_TEST( qoiValue );
    CPPUNIT_TEST_SUITE_END();

    public:
      Communicator comm, physicsComm;
      GetPot inputfile;
      std::shared_ptr<PhysicsLibmesh<T_S,T_P> > myPhysics ; 
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


        inputfile.set("solutions","qoi") ;
        inputfile.set("nElem",1) ;

        myPhysics = std::shared_ptr<PhysicsLibmesh<T_S,T_P> > (
            new PhysicsViscousBurgers<T_S,T_P>(physicsComm,inputfile ) 
            );

        dimension = 2;
        myParameters.reserve(dimension);
        myParameters.push_back( 
            std::shared_ptr<AGNOS::Parameter>(
              new AGNOS::Parameter(UNIFORM,-1.0, 1.0) )
          ); 
        myParameters.push_back( 
            std::shared_ptr<AGNOS::Parameter>(
              new AGNOS::Parameter(UNIFORM,-1.0, 1.0) )
          );
      }

      void tearDown( )
      {
      }

      
      void qoiValue()
      {
        std::cout 
          << "\n--------------------------\n "
          << "Testing ViscousBurgers for evaluation of QoI"
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
        paramValue.zero();

        T_P testValue; 

        mySurrogate->build( );
        testValue = mySurrogate->evaluate( "qoi", paramValue );

        std::cout << "testValue = " << testValue(0) << std::endl;
        std::cout << " n_elem = " << myPhysics->getMesh( ).n_active_elem() << std::endl;

        CPPUNIT_ASSERT( 
            std::abs(testValue(0) - 10.0) <= 1e-9
            );

        delete mySurrogate;
      }



  };

  CPPUNIT_TEST_SUITE_REGISTRATION( ViscousBurgersTest );

//________________________________________________________________//

// EOF
