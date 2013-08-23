

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Parallel
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// local includes
#include <iostream>
#include <stdio.h>
#include <assert.h>

#include "agnosDefines.h"
#include <mpi.h>
#include "Parameter.h"
#include "PseudoSpectralTensorProduct.h"
#include "Parameter.h"

#include "PhysicsCatenary.h"

using namespace AGNOS;

//________________________________________________________________//

typedef libMesh::DenseVector<double> T_P ;
typedef libMesh::DenseVector<double> T_S ;


int already_initialized;
int already_finalized;
Communicator comm;
GetPot inputfile = GetPot();

void finalize_mpi(void)
{
  MPI_Finalized(&already_finalized);
  if(!already_finalized)
    MPI_Finalize();
}

void initialize_mpi()
{
  MPI_Initialized(&already_initialized);
  if (!already_initialized)
  {
    MPI_Init(NULL,NULL);
    MPI_Initialized(&already_initialized);
    atexit(finalize_mpi);
  }
}


BOOST_AUTO_TEST_SUITE(Catenary_tensorProduct)


  unsigned int dimension = 1;
  
  std::vector<Parameter*> myParameters(
      dimension, 
      new Parameter(UNIFORM, 1.0,3.0)
      ); 

  PhysicsCatenary<T_S,T_P>* myPhysics = new PhysicsCatenary<T_S,T_P>(
      comm, inputfile );

BOOST_AUTO_TEST_CASE(Catenary_N0)
{

  initialize_mpi();
  comm = MPI_COMM_WORLD;
  /* inputfile.set("physics/solutions","primal") ; */

  if (comm.rank()==0)
    BOOST_TEST_MESSAGE(" Testing Catenary with N=0");

  std::vector<unsigned int> myOrder(dimension,0);
  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
    new PseudoSpectralTensorProduct<T_S,T_P>(
        comm,
        myPhysics, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::map< std::string, std::vector<T_P> > myCoeff ;

  if ( 0 % comm.size() == comm.rank() )
  {
    myCoeff = mySurrogate->getCoefficients();
    BOOST_CHECK_CLOSE( myCoeff["primal"][0](0) , -10.0/16.0, 1e-9 );
  }

}

BOOST_AUTO_TEST_CASE(Catenary_N1)
{
  initialize_mpi();
  comm = MPI_COMM_WORLD;
  /* inputfile.set("physics/solutions","primal") ; */

  if (comm.rank()==0)
    BOOST_TEST_MESSAGE(" Testing Catenary with N=1");

  std::vector<unsigned int> myOrder(dimension,1);
  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
    new PseudoSpectralTensorProduct<T_S,T_P>(
        comm,
        myPhysics, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::map< std::string, std::vector<T_P> > myCoeff;

  myCoeff = mySurrogate->getCoefficients();

  if (comm.rank() == 0)
  {
    // 12/11 is the sum of u(\xi_j) 
    BOOST_CHECK_CLOSE( 
        myCoeff["primal"][0](0),
        -10.0/16.0 * (12./11.), 1e-9 );

    // -2 sqrt(3)/11 is correct contribution from poly evals
    BOOST_CHECK_CLOSE( 
        myCoeff["primal"][1](0) , 
        -10.0/16.0 * ( -2.0 * std::sqrt(3.0) / 11.0 ), 
        1e-9 );
  }

}

BOOST_AUTO_TEST_CASE(Catenary_N4)
{
  initialize_mpi();
  comm = MPI_COMM_WORLD;
  /* inputfile.set("physics/solutions","primal") ; */

  if (comm.rank()==0)
    BOOST_TEST_MESSAGE(" Testing Catenary with N=4");

  std::vector<unsigned int> myOrder(dimension,4);
  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
    new PseudoSpectralTensorProduct<T_S,T_P>(
        comm,
        myPhysics, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::map< std::string, std::vector<T_P> > myCoeff ;

  myCoeff = mySurrogate->getCoefficients();

  if (comm.rank() == 0)
  {
    // Coefficients generated from pmpack for comparison
    BOOST_CHECK_CLOSE( 
        myCoeff["primal"][0](0) , 
        -6.866307761327950e-01,
        1e-9 );
    std::cout << "myCoeff[0]=" << myCoeff["primal"][0](0) << std::endl;
    
    BOOST_CHECK_CLOSE( 
        myCoeff["primal"][1](0) , 
        2.134952711438086e-01,
        1e-9 );
    std::cout << "myCoeff[1]=" << myCoeff["primal"][1](0) << std::endl;
   
    BOOST_CHECK_CLOSE( 
        myCoeff["primal"][2](0) , 
        -5.918708419582225e-02, 
        1e-9 );
    std::cout << "myCoeff[2]=" << myCoeff["primal"][2](0) << std::endl;

    BOOST_CHECK_CLOSE( 
        myCoeff["primal"][3](0) , 
        1.602406581398491e-02,
        1e-9 );
    std::cout << "myCoeff[3]=" << myCoeff["primal"][3](0) << std::endl;

    BOOST_CHECK_CLOSE( 
        myCoeff["primal"][4](0) , 
        -4.037685060565489e-03,
        1e-9 );
    std::cout << "myCoeff[4]=" << myCoeff["primal"][4](0) << std::endl;
  }

}

BOOST_AUTO_TEST_CASE(Catenary_mean)
{
  initialize_mpi();
  comm = MPI_COMM_WORLD;
  /* inputfile.set("physics/solutions","primal") ; */

  if (comm.rank()==0)
    BOOST_TEST_MESSAGE(" Testing Catenary with N=4 for mean calculation");

  std::vector<unsigned int> myOrder(dimension,4);
  PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = 
    new PseudoSpectralTensorProduct<T_S,T_P>(
        comm,
        myPhysics, 
        myParameters, 
        myOrder  
        );

  mySurrogate->build( );
  std::map< std::string, T_P > myMean ;

  myMean = mySurrogate->mean();

  // Coefficients generated from pmpack for comparison
  BOOST_CHECK_CLOSE( 
      myMean["primal"](0) , 
      -6.866307761327950e-01,
      1e-9 );
    

}

BOOST_AUTO_TEST_CASE(Catenary_convergence)
{
  initialize_mpi();
  comm = MPI_COMM_WORLD;
  /* inputfile.set("physics/solutions","primal") ; */

  if (comm.rank()==0)
    BOOST_TEST_MESSAGE(" Testing Catenary convergence");

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
    testValue = mySurrogate->evaluate( "primal", paramValue ) ;
  }

  if (comm.rank() == 0)
    BOOST_CHECK_CLOSE(  testValue(0) , -10.0/(8.0 * paramValue(0) ), 1e-9 );


}


BOOST_AUTO_TEST_SUITE_END()


//________________________________________________________________//

// EOF
