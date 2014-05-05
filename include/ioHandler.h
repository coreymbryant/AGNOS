
#ifndef IO_HANDLER_H
#define IO_HANDLER_H

#include "Driver.h"
#include <H5Cpp.h>

using namespace H5;

namespace AGNOS{

  // TODO USE user-block to tag date, time, commit ref etc
  // defined in file creation/access property list
  // H5P_DEFAULT used for now

  /** write parameter info: dim, type, bounds; to HDF5 file*/
  // H5F_ACC_TRUNC
  void writeParameters(
      std::string fileName,
      std::vector<std::shared_ptr<AGNOS::Parameter> > parameters);

  /** read parameter info: dim, type, bounds; from HDF5 file*/
  // H5F_RDONLY
  void readParameters(
      std::string fileName,
      std::vector<std::shared_ptr<AGNOS::Parameter> > parameters);

  /** write surrogate properties to HDF5 file */
  // H5F_ACC_TRUNC
  void writeSurrogate( 
      std::string fileName,
      std::shared_ptr<SurrogateModel<T_S,T_P> > surrogate );

  /** read surrogate properties from HDF5 file */
  // H5F_RDONLY
  void readSurrogate( 
      std::string fileName,
      std::shared_ptr<SurrogateModel<T_S,T_P> > surrogate );

  // TODO: do we need physics read and write or do we rely on input file
  //  - would at least need mesh info to make sure it matches dimension of coeff
  //  vectors

}
#endif // IO_HANDLER_H
