
#ifndef IO_HANDLER_H
#define IO_HANDLER_H

#include "Driver.h"
#include <H5Cpp.h>

using namespace H5;

namespace AGNOS{


  // TODO USE user-block to tag date, time, commit ref etc
  // defined in file creation/access property list
  // H5P_DEFAULT used for now
  
  /** write user block to file for reference purposes */
  H5File* writeUserBlock( std::string fileName );

  /** read in user block and get a reference to H5File object */
  // H5F_RDONLY
  H5File* readUserBlock( std::string fileName );

  /** write parameter info: dim, type, bounds; to HDF5 file*/
  void writeParameters(
      CommonFG* common,
      std::vector<std::shared_ptr<AGNOS::Parameter> > parameters);

  /** read parameter info: dim, type, bounds; from HDF5 file*/
  // H5F_RDONLY
  void readParameters(
      CommonFG* common,
      std::vector<std::shared_ptr<AGNOS::Parameter> >& parameters);

  /** write surrogate properties to HDF5 file */
  void writeSurrogate( 
      CommonFG* common,
      std::shared_ptr<SurrogateModel<T_S,T_P> >& surrogate );

  /** read surrogate properties from HDF5 file */
  // H5F_RDONLY
  void readSurrogate( 
      CommonFG* common,
      std::vector<unsigned int>& order,
      std::set<std::string>& computeSolutions,
      std::vector< std::vector<unsigned int> >& indexSet,
      std::map< std::string, LocalMatrix >& coefficients
      );

  // TODO: do we need physics read and write or do we rely on input file
  //  - would at least need mesh info to make sure it matches dimension of coeff
  //  vectors

}
#endif // IO_HANDLER_H
