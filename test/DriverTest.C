// hack to allow access to private/protected members
#define private public
#define protected public

// local
#include "agnosDefines.h"
#include "Driver.h"

// boost
#define BOOST_TEST_MODULE driverTest
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace AGNOS;

BOOST_AUTO_TEST_CASE( driver_constructor )
{
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];
  char* name = new char[27];
  strcpy(name, "test");
  av[0] = name;
  MPI_Init(&ac,&av);
  Communicator comm(MPI_COMM_WORLD);
  LibMeshInit libmesh_init(ac, av, comm.get()) ;
  delete av, name;

  GetPot input( "driver.in" );
  AGNOS::Driver driver( comm, comm, input );

  BOOST_REQUIRE( (driver._maxIter == 1) );
  BOOST_REQUIRE( (driver._adaptiveDriver == false) ) ;
  BOOST_REQUIRE( (driver._refinePercentage == 0.20) );
  BOOST_REQUIRE( (driver._simultRefine == false) );
  BOOST_REQUIRE( (driver._paramDim == 1) );
  BOOST_REQUIRE( (driver._nInitialHRefinements == 0) );

  BOOST_REQUIRE( (driver._refinePhysics == true) );
  BOOST_REQUIRE( (driver._uniformRefine == true) );

  BOOST_REQUIRE( (driver._surrogateNames.size() == 1) );
  BOOST_REQUIRE( 
      (driver._surrogateNames.begin()->first == "primarySurrogate") );
  BOOST_REQUIRE( 
      (driver._surrogateNames.begin()->second == 0 ) );
  BOOST_REQUIRE( (driver._refineSurrogate == false) );
  BOOST_REQUIRE( (driver._hRefine == false) );
  BOOST_REQUIRE( (driver._pRefine == true) );
  BOOST_REQUIRE( (driver._anisotropic == false) );
  BOOST_REQUIRE( (driver._pIncrement.size() == 1) );
  BOOST_REQUIRE( (driver._pIncrement[0] == 1) );

  BOOST_REQUIRE( (driver._elemsToUpdate.size() == 1) );

  BOOST_REQUIRE( (driver._outputFilename == "./testFile.dat") );
  BOOST_REQUIRE( (driver._solutionsToPrint.size() == 3) );
  BOOST_REQUIRE( (driver._solutionsToPrint[0] == "primal") );
  BOOST_REQUIRE( (driver._solutionsToPrint[1] == "adjoint") );
  BOOST_REQUIRE( (driver._solutionsToPrint[2] == "qoi") );
  BOOST_REQUIRE( (driver._computeMeans        == true) );
  BOOST_REQUIRE( (driver._computeNorms        == true) );
  BOOST_REQUIRE( (driver._outputIterations    == false) );
  BOOST_REQUIRE( (driver._outputCoefficients  == true) );
  BOOST_REQUIRE( (driver._outputWeights       == true) );
  BOOST_REQUIRE( (driver._outputPoints        == true) );
  BOOST_REQUIRE( (driver._outputIndexSet      == true) );

  BOOST_REQUIRE( (driver._generateSamples     == false) );
  BOOST_REQUIRE( (driver._sampleFile          == "./sampleFile") );
  BOOST_REQUIRE( (driver._nSamples            == 10000) );

}

BOOST_AUTO_TEST_CASE( driver_build )
{
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];
  char* name = new char[27];
  strcpy(name, "test");
  av[0] = name;
  MPI_Init(&ac,&av);
  Communicator comm(MPI_COMM_WORLD);
  LibMeshInit libmesh_init(ac, av, comm.get()) ;
  delete av, name;

  GetPot input( "driver.in" );
  AGNOS::Driver driver( comm, comm, input );

  driver.run();
  
  BOOST_REQUIRE( (driver._elemsToUpdate.empty()) );

}

BOOST_AUTO_TEST_CASE( driver_build_evaluator )
{
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];
  char* name = new char[27];
  strcpy(name, "test");
  av[0] = name;
  MPI_Init(&ac,&av);
  Communicator comm(MPI_COMM_WORLD);
  LibMeshInit libmesh_init(ac, av, comm.get()) ;
  delete av, name;

  GetPot input ;
  input = GetPot( );
  input.set("driver/evaluator",true);
  input.set("driver/evaluatorFile","test.h5");
  AGNOS::Driver driver( comm, comm, input );

  
  BOOST_REQUIRE( (driver._activeElems.size() == 16) );

}

BOOST_AUTO_TEST_CASE( driver_build_from_restart )
{
}
