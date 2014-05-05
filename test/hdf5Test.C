

// local
#include "ioHandler.h"

// boost
#define BOOST_TEST_MODULE hdf5Test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace AGNOS;

BOOST_AUTO_TEST_CASE( ioHandler_user_block )
{
  std::string fileName = "test.h5";
  writeUserBlock( fileName );
}

BOOST_AUTO_TEST_CASE( ioHadnler_write_parameters )
{
  std::string fileName = "test.h5";
  writeUserBlock( fileName );

  // construct parameter object
  unsigned int dimension = 3;
  std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;
  for(unsigned int p=0;p<dimension;p++)
    parameters.push_back(
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter("CONSTANT",1.0,1.0) )
        );


  // write parameters to file
  writeParameters( fileName, parameters );

  //read in parameters from file
  std::vector<std::shared_ptr<AGNOS::Parameter> > newParameters;
  readParameters( fileName, newParameters );


  BOOST_REQUIRE( (parameters.size() == newParameters.size() ) );
  for(unsigned int i = 0; i<parameters.size();i++)
  {
    BOOST_REQUIRE( (parameters[i]->type() == newParameters[i]->type()) );
    BOOST_REQUIRE( (parameters[i]->min() == newParameters[i]->min()) );
    BOOST_REQUIRE( (parameters[i]->max() == newParameters[i]->max()) );
  }

}
