

// local
#include "ioHandler.h"

// boost
#define BOOST_TEST_MODULE hdf5Test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace AGNOS;

BOOST_AUTO_TEST_CASE( hdf5_open_close )
{
  std::string fileName = "test.h5";

  try {
    // create a new file using default properties
    H5File h5_file(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    // close file
    h5_file.close();
  }

  catch( FileIException error )
  {
    error.printError();
    exit(1);
  }


}

BOOST_AUTO_TEST_CASE( ioHadnler_write_parameters )
{
  std::string fileName = "test.h5";

  // construct parameter object
  unsigned int dimension = 3;
  std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;
  for(unsigned int p=0;p<dimension;p++)
    parameters.push_back(
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter("CONSTANT",1.0,1.0) )
        );

  writeParameters( fileName, parameters );



}
