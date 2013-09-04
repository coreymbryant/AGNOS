

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

// local includes
#include <iostream>
#include <stdio.h>
#include <assert.h>

#include "agnosDefines.h"
#undef AGNOS_DEBUG
#define AGNOS_DEBUG 1

#include "Parameter.h"
#include "PseudoSpectralTensorProduct.h"

using namespace AGNOS;

//________________________________________________________________//

  typedef libMesh::DenseVector<double> T_P ;
  typedef libMesh::DenseVector<double> T_S ;
  // linear test function
  void linearFunction (
      const T_S& paramVec, 
      std::map<std::string,T_P>& solutionVectors
      )
  {
    T_P returnVec(1);
    returnVec(0)=paramVec(0);
    for(unsigned int dim=1; dim<paramVec.size(); dim++)
      returnVec(0) *= paramVec(dim) ;

    solutionVectors.clear();
    solutionVectors.insert(std::pair<std::string,T_P>("primal",returnVec) );
  }


  // quadratic test function
  void quadFunction (
      const T_S& paramVec, 
      std::map<std::string,T_P>& solutionVectors
      )
  {
    T_P returnVec(1);
    returnVec(0) = paramVec(0) * paramVec(0) ;
    for(unsigned int dim=1; dim<paramVec.size(); dim++)
      returnVec(0) *= (paramVec(dim)*paramVec(dim)) ;

    solutionVectors.clear();
    solutionVectors.insert(std::pair<std::string,T_P>("primal",returnVec) );
  }


  // mixed-order 5 dimensional exammple
  void mixedFunction (
      const T_S& paramVec, 
      std::map<std::string,T_P>& solutionVectors
      )
  {
    T_P returnVec(1);
    returnVec(0) = paramVec(2)  + paramVec(1)*paramVec(1)*paramVec(0) + std::pow(paramVec(3),3) ;

    solutionVectors.clear();
    solutionVectors.insert(std::pair<std::string,T_P>("primal",returnVec) );
  }


  class PolynomialTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( PolynomialTest );
    CPPUNIT_TEST( Linear_1D );
    CPPUNIT_TEST( Linear_ND );
    CPPUNIT_TEST( Quad_1D );
    CPPUNIT_TEST( Quad_ND );
    CPPUNIT_TEST( OrderN_1D );
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


      void Linear_1D()
      {
        std::cout << " Testing 1D linear approximation" << std::endl;

        unsigned int dimension = 1;
        std::vector<unsigned int> myOrder(dimension,3);

        std::vector<std::shared_ptr<AGNOS::Parameter> > myParameters(
            dimension, 
            std::shared_ptr<AGNOS::Parameter>(
              new AGNOS::Parameter(UNIFORM, -1.0,1.0)
              )
            ); 

        PhysicsModel<T_S,T_P>* myPhysics=
          new PhysicsModel<T_S,T_P>(physicsComm,inputfile) ;
        myPhysics->attach_compute_function(&linearFunction);

        PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
          PseudoSpectralTensorProduct<T_S,T_P>(
              comm,
              myPhysics, 
              myParameters, 
              myOrder  
              );

        mySurrogate->build( );
        std::map< std::string, LocalMatrix > myCoeff
          = mySurrogate->getCoefficients( );

        if (comm.rank() == 0)
          CPPUNIT_ASSERT( 
            std::abs( myCoeff["primal"](1,0) - std::sqrt(1.0/3.0) ) <= 1e-9 
          );

        delete mySurrogate;
        delete myPhysics;

      }


      void Linear_ND()
      {
        std::cout << " Testing 7D linear approximation" << std::endl;

        unsigned int dimension = 7;
        std::vector<unsigned int> myOrder(dimension,1);
        std::vector<std::shared_ptr<AGNOS::Parameter> > myParameters(
            dimension, 
            std::shared_ptr<AGNOS::Parameter>(
              new AGNOS::Parameter(UNIFORM, -1.0,1.0)
              )
            ); 

        PhysicsModel<T_S,T_P>* myPhysics=
          new PhysicsModel<T_S,T_P>(physicsComm,inputfile) ;
        myPhysics->attach_compute_function(&linearFunction);

        PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
          PseudoSpectralTensorProduct<T_S,T_P>(
              comm,
              myPhysics, 
              myParameters, 
              myOrder  
              );

        mySurrogate->build( );
        std::map< std::string, LocalMatrix > myCoeff
          = mySurrogate->getCoefficients( );


        if (comm.rank() == 0)
          CPPUNIT_ASSERT( 
              std::abs( 
                myCoeff["primal"](myCoeff["primal"].m()-1,0) 
                - 
                std::pow(static_cast<double>(std::sqrt(1.0/3.0)),
                  static_cast<int>(dimension))
                ) <= 1e-4 );

        delete mySurrogate;
        delete myPhysics;
      }

      void Quad_1D()
      {
        std::cout << " Testing 1D quadratic approximation" << std::endl;

        unsigned int dimension = 1;
        std::vector<unsigned int> myOrder(dimension,2);

        std::vector<std::shared_ptr<AGNOS::Parameter> > myParameters(
            dimension, 
            std::shared_ptr<AGNOS::Parameter>(
              new AGNOS::Parameter(UNIFORM, -1.0,1.0)
              )
            ); 

        PhysicsModel<T_S,T_P>* myPhysics=
          new PhysicsModel<T_S,T_P>(physicsComm,inputfile) ;
        myPhysics->attach_compute_function(&quadFunction);

        PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
          PseudoSpectralTensorProduct<T_S,T_P>(
              comm,
              myPhysics, 
              myParameters, 
              myOrder  
              );

        mySurrogate->build( );
        std::map< std::string, LocalMatrix > myCoeff
          = mySurrogate->getCoefficients( );

        if (comm.rank() == 0)
        {
          CPPUNIT_ASSERT( std::abs( myCoeff["primal"](0,0) - 1.0/3.0 ) <=  1e-9 );
          CPPUNIT_ASSERT( 
              std::abs( myCoeff["primal"](2,0) - 2.0/3.0 * std::sqrt(1.0/5.0) ) 
              <= 1e-9 );
        }

        delete mySurrogate;
        delete myPhysics;
      }

      void Quad_ND()
      {
        std::cout << " Testing 3D quadratic approximation" << std::endl;

        unsigned int dimension = 3;
        std::vector<unsigned int> myOrder(dimension,2);

        std::vector<std::shared_ptr<AGNOS::Parameter> > myParameters(
            dimension, 
            std::shared_ptr<AGNOS::Parameter>(
              new AGNOS::Parameter(UNIFORM, -1.0,1.0)
              )
            ); 

        PhysicsModel<T_S,T_P>* myPhysics=
          new PhysicsModel<T_S,T_P>(physicsComm,inputfile) ;
        myPhysics->attach_compute_function(&quadFunction);

        PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
          PseudoSpectralTensorProduct<T_S,T_P>(
              comm,
              myPhysics, 
              myParameters, 
              myOrder  
              );

        mySurrogate->build( );
        std::map< std::string, LocalMatrix > myCoeff
          = mySurrogate->getCoefficients( );

        if (comm.rank() == 0)
        {
          CPPUNIT_ASSERT( 
              std::abs( 
                myCoeff["primal"](0,0)  - 
                std::pow(static_cast<double>(1.0/3.0), static_cast<int>(dimension)) 
                ) <= 1e-9 );
          CPPUNIT_ASSERT( 
              std::abs( myCoeff["primal"](myCoeff["primal"].m()-1,0) - 
                std::pow( static_cast<double>(std::sqrt(1.0/5.0) * 2.0/3.0),
                static_cast<int>(dimension)) 
                ) <= 1e-9 );
        }
        delete mySurrogate;
        delete myPhysics;
      }

      void OrderN_1D()
      {
        std::cout << " Testing mixed-order 5 dimensional example" << std::endl;

        unsigned int dimension = 5;
        std::vector<unsigned int> myOrder(dimension,0);
        myOrder[0] = 1;
        myOrder[1] = 2;
        myOrder[2] = 1;
        myOrder[3] = 3;
        myOrder[4] = 0;

        std::vector<std::shared_ptr<AGNOS::Parameter> > myParameters(
            dimension, 
            std::shared_ptr<AGNOS::Parameter>(
              new AGNOS::Parameter(UNIFORM, -1.0,1.0)
              )
            ); 

        PhysicsModel<T_S,T_P>* myPhysics=
          new PhysicsModel<T_S,T_P>(physicsComm,inputfile) ;
        myPhysics->attach_compute_function(&mixedFunction);

        PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
          PseudoSpectralTensorProduct<T_S,T_P>(
              comm,
              myPhysics, 
              myParameters, 
              myOrder  
              );

        mySurrogate->build( );
        std::map< std::string, LocalMatrix > myCoeff
          = mySurrogate->getCoefficients( );

        /* mySurrogate->printIntegrationWeights() ; */
        /* mySurrogate->printIntegrationPoints() ; */
        /* mySurrogate->printIndexSet() ; */
        /* for(unsigned int coeff=0; coeff<myCoeff.size(); coeff++) */
        /*   for(unsigned int dim=0; dim<myCoeff[coeff].size(); dim++) */
        /*     std::cout << "coeff[" << coeff << "](" << dim << ") = " */ 
        /*       << myCoeff[coeff](dim) << std::endl; */

        if (comm.rank() == 0)
        {
          // 0 0 1 0 0  = 1 * 1/norm(psi) from normailziation
          CPPUNIT_ASSERT( 
              std::abs( myCoeff["primal"](4,0)  - std::sqrt(1.0/3.0) )
              <= 1e-9 );

          // 1 2 0 0 0  = 1
          // have to use combination of terms since psi_2 = 1/2*(3x^2-1) 
          CPPUNIT_ASSERT(
              std::abs( myCoeff["primal"](40,0) - 2.0/(3.0 * std::sqrt(3.0*5.0)) )
              <= 1e-9 );
          CPPUNIT_ASSERT(
              std::abs( myCoeff["primal"](24,0) - 1.0 / ( 3.0 * std::sqrt(3.0) ) )
              <= 1e-9 );

          // 0 0 0 3 0  = 1
          // have to use combination of terms since psi_2 = 1/2*(5x^3-3x) 
          CPPUNIT_ASSERT(
              std::abs( myCoeff["primal"](3,0) - 2.0 / ( 5.0 * std::sqrt(7.0) ) )
              <= 1e-9 );
          CPPUNIT_ASSERT(
               std::abs( myCoeff["primal"](1,0) - 3.0 / ( 5.0 * std::sqrt(3.0) ))
               <= 1e-9 );
        }
        delete mySurrogate;
        delete myPhysics;
      }


  };

  CPPUNIT_TEST_SUITE_REGISTRATION( PolynomialTest );
//________________________________________________________________//

// EOF
