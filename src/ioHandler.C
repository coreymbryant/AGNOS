
#include "ioHandler.h"

namespace AGNOS{

  /** write parameter info: dim, type, bounds; to HDF5 file*/
  // H5F_ACC_TRUNC
  void writeParameters(
      std::string fileName,
      std::vector<std::shared_ptr<AGNOS::Parameter> > parameters)
  {

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

      return;
    }

  /** read parameter info: dim, type, bounds; from HDF5 file*/
  // H5F_RDONLY
  void readParameters(
      std::string fileName,
      std::vector<std::shared_ptr<AGNOS::Parameter> > parameters)
  {

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

    return;
  }

  /** write surrogate properties to HDF5 file */
  // H5F_ACC_TRUNC
  void writeSurrogate( 
      std::string fileName,
      std::shared_ptr<SurrogateModel<T_S,T_P> > surrogate )
  {

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

    return;
  }

  /** read surrogate properties from HDF5 file */
  // H5F_RDONLY
  void readSurrogate( 
      std::string fileName,
      std::shared_ptr<SurrogateModel<T_S,T_P> > surrogate )
  {

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

    return;
  }


} // end namespace
