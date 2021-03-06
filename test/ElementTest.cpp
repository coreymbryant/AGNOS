

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

  // constant test function
  void constantFunction (
      std::set<std::string>& computeSolutions,
      const T_S& paramVec, 
      std::map<std::string,T_P>& solutionVectors
      )
  {
    T_P returnVec(1);
    returnVec(0)= 2.0;
    solutionVectors.clear();
    solutionVectors.insert(std::pair<std::string,T_P>("primal",returnVec) );
  }



// Element test class
  class ElementTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( ElementTest );
    CPPUNIT_TEST( coversTest  );
    CPPUNIT_TEST( splitSamePoly );
    CPPUNIT_TEST( splitConstant );
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

      /** Split element and test which element contains new eval point  */
      void coversTest()
      {
        std::cout 
          << "\n--------------------------\n "
          << "Testing Element framework for covering eval point "
          << std::endl;

        myParameters.push_back( 
            std::shared_ptr<AGNOS::Parameter>(
              new AGNOS::Parameter(CONSTANT,1.0, 1.0) )
            );
        dimension++;

        std::vector<unsigned int> myOrder(dimension,1);

        std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > > surrogates ;
        surrogates.push_back( 
            std::shared_ptr<SurrogateModelBase<T_S,T_P> > (
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

        T_S paramValue0(dimension);
        T_S paramValue1(dimension);
        T_S paramValue2(dimension);
        paramValue0(0) = 1.25;
        paramValue0(1) = 1.0;
        paramValue1(0) = 2.25;
        paramValue1(1) = 1.0;
        paramValue2(0) = 2.25;
        paramValue2(1) = 2.0;


        std::vector< Element<T_S,T_P> > children = baseElement.split() ;
        for (unsigned int c=0; c<children.size(); c++)
        {
          bool test0 = children[c].covers( paramValue0 );
          bool test1 = children[c].covers( paramValue1 );
          std::cout << "child: " << c << " min = " <<
            children[c].parameters()[0]->min() << std::endl;
          std::cout << "child: " << c << " max = " <<
            children[c].parameters()[0]->max() << std::endl;
          std::cout << "child: " << c << " covers0? = " << test0 << std::endl;
          std::cout << "child: " << c << " covers1? = " << test1 << std::endl;

          if ( c == 0)
            CPPUNIT_ASSERT( test0 );
          if ( c == 1)
            CPPUNIT_ASSERT( test1 );

        }

        bool test3 = children[0].covers( paramValue2 );
        CPPUNIT_ASSERT( !test3 );



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

        std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > > surrogates ;
        surrogates.push_back( 
            std::shared_ptr<SurrogateModelBase<T_S,T_P> > (
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
              std::abs(testValue(0) - baseValue(0)) <= 1e-12
              );

        }




      }

      /********************************************//**
       * \brief Test that splitting a constant valued function doesn't actually
       * change anything, i.e. each *new* surrogate still computes correct value
       *
       * 
       ***********************************************/
      void splitConstant()
      {
        std::cout 
          << "\n--------------------------\n "
          << "Testing Element framework when splitting: constant "
          << std::endl;

        std::vector<unsigned int> myOrder(dimension,2);

        std::shared_ptr<PhysicsModel<T_S,T_P> > myPhysics;
        std::shared_ptr<PhysicsUser<T_S,T_P> > constantPhysics( 
            new PhysicsUser<T_S,T_P>(physicsComm,inputfile) 
            );
        constantPhysics->attach_compute_function(&constantFunction);
        myPhysics = constantPhysics;

        std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > > surrogates ;
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
        T_P testValue, testMean ;

        T_P baseValue = baseElement.surrogates()[0]->evaluate( "primal", paramValue );
        std::cout << "base testValue = " << baseValue(0) << std::endl;

        std::map<std::string,T_P> baseMeans =
          baseElement.surrogates()[0]->mean( );
        T_P baseMean( baseMeans["primal"] );
        std::cout << "base meanValue = " << baseMean(0) << std::endl;

        double baseWeight = baseElement.weight();
        std::cout << "base weight = " << baseWeight << std::endl;
        double testWeight = 0.;

        std::vector< Element<T_S,T_P> > children = baseElement.split() ;
        for (unsigned int c=0; c<children.size(); c++)
        {
          paramValue(0) = 0.5 * ( children[c].parameters()[0]->min() 
            + children[c].parameters()[0]->max() );
          std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > >
            childSurrogates;
          childSurrogates.push_back( 
              std::shared_ptr<SurrogateModel<T_S,T_P> > (
                new PseudoSpectralTensorProduct<T_S,T_P>(
                    comm,
                    children[c].physics(),
                    children[c].parameters(), 
                    myOrder  
                    )
                )
            );
          children[c].setSurrogates( childSurrogates );
          children[c].surrogates()[0]->build();

          testValue = children[c].surrogates()[0]->evaluate( "primal", paramValue );
          std::cout << "child: " << c << " testValue = " << testValue(0) << std::endl;
          std::cout << "child: " << c << " diff = " << 
              std::abs(testValue(0) - baseValue(0)) << std::endl;
          std::map<std::string,T_P> childMeans =
            children[c].surrogates()[0]->mean( );
          testMean = childMeans["primal"] ;
          std::cout << "child: " << c << " meanValue = " << testMean(0) <<
            std::endl;

          std::cout << "child: " << c << " weight = " << children[c].weight() <<
            std::endl;
          testWeight += children[c].weight();

          CPPUNIT_ASSERT( 
              std::abs(testValue(0) - baseValue(0)) <= 1e-12
              );
          CPPUNIT_ASSERT( 
              std::abs(testMean(0) - baseMean(0)) <= 1e-12
              );

        }

        std::cout << " sum of child weights = " << testWeight << std::endl;
        CPPUNIT_ASSERT( 
            std::abs(testWeight - baseWeight) <= 1e-12
            );




      }



  };

  CPPUNIT_TEST_SUITE_REGISTRATION( ElementTest );

//________________________________________________________________//

// EOF
