

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
#include "Element.h"

#include "PhysicsCatenary.h"


using namespace AGNOS;

//________________________________________________________________//

  typedef libMesh::DenseVector<double> T_P ;
  typedef libMesh::DenseVector<double> T_S ;
  // linear test function

  class ElementTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( ElementTest );
    CPPUNIT_TEST( splitSamePoly );
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


        myPhysics = std::shared_ptr<PhysicsModel<T_S,T_P> > (
            new PhysicsCatenary<T_S,T_P>(physicsComm,inputfile ) 
            );

        dimension = 1;
        myParameters.reserve(dimension);
        myParameters.push_back( 
            std::shared_ptr<AGNOS::Parameter>(
              new AGNOS::Parameter(UNIFORM,1.0, 3.0) )
          ); 
      }

      void tearDown( )
      {
      }

      
      /********************************************//**
       * \brief Test to confirm that splitting an Element, but not rebuilding
       * surrogate maintains the same surrogate model. This is used when we want
       * to refine one child but keep the parent surrogate model for other
       * children to save on computations
       *
       * 
       ***********************************************/
      void splitSamePoly()
      {
        std::cout 
          << "\n--------------------------\n "
          << "Testing Element framework when splitting "
          << std::endl;

        std::vector<unsigned int> myOrder(dimension,3);

        std::vector< std::shared_ptr<SurrogateModel<T_S,T_P> > > surrogates ;
        surrogates.push_back( 
            std::shared_ptr<SurrogateModel<T_S,T_P> > (
              new PseudoSpectralTensorProduct<T_S,T_P>(
                  comm,
                  myPhysics,
                  myParameters, 
                  myOrder  
                  )
              )
          );

        AGNOS::Element<T_S,T_P> baseElement(
            myParameters,
            surrogates,
            myPhysics
            );
        baseElement.surrogates()[0]->build( );

        T_S paramValue(dimension);
        paramValue(0) = 1.5;
        T_P testValue ;

        T_P baseValue = baseElement.surrogates()[0]->evaluate( "primal", paramValue );
        std::cout << "base testValue = " << baseValue(0) << std::endl;

        std::vector< Element<T_S,T_P> > children = baseElement.split() ;
        for (unsigned int c=0; c<children.size(); c++)
        {
          testValue = children[c].surrogates()[0]->evaluate( "primal", paramValue );
          std::cout << "child: " << c << " testValue = " << testValue(0) << std::endl;

          CPPUNIT_ASSERT( 
              std::abs(testValue(0) - baseValue(0)) <= 1e-16
              );

        }




      }



  };

  CPPUNIT_TEST_SUITE_REGISTRATION( ElementTest );

//________________________________________________________________//

// EOF
