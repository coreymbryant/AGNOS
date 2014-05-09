
#ifndef IO_HANDLER_H
#define IO_HANDLER_H

#include "Driver.h"
#include <H5Cpp.h>

using namespace H5;

namespace AGNOS{

  class H5IO 
  {
    public:
      /**  constructor */
      H5IO( std::string fileName, int accessFlag = H5F_ACC_EXCL );
      /** Default destructor */
      ~H5IO();

      void close();

      // TODO USE user-block to tag date, time, commit ref etc
      // defined in file creation/access property list
      // H5P_DEFAULT used for now
      /** write user block to file for reference purposes */
      void writeUserBlock(  );

      /** read in user block and get a reference to H5File object */
      // H5F_RDONLY
      void readUserBlock(  );

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
          std::shared_ptr<SurrogateModelBase<T_S,T_P> >& surrogate,
          std::string surrName = "surrogate" );

      /** read surrogate properties from HDF5 file */
      // H5F_RDONLY
      void readSurrogate( 
          CommonFG* common,
          std::vector<unsigned int>& order,
          std::set<std::string>& computeSolutions,
          std::vector< std::vector<unsigned int> >& indexSet,
          std::map< std::string, LocalMatrix >& coefficients,
          std::string surrName = "surrogate"
          );

      // TODO: do we need physics read and write or do we rely on input file
      //  - would at least need mesh info to make sure it matches dimension of
      //  coeff vectors
      /** write physics data */
      void writePhysics(
          CommonFG* common,
          std::shared_ptr<PhysicsModel<T_S,T_P> >&  physics );
      /** read physics data */
      void readPhysics(
          CommonFG* common,
          std::shared_ptr<PhysicsModel<T_S,T_P> >&  physics );


      /** write element */
      void writeElement(
          CommonFG* common,
          AGNOS::Element<T_S,T_P>& element ,
          unsigned int elemId = 0);
      /** read element */
      void readElement(
          CommonFG* common,
          AGNOS::Element<T_S,T_P>& element,
          unsigned int elemId = 0);


      /** write simulation data */
      /** read simulation data */

      /** reference to underlying H5File object */
      H5File* file() const { return _h5File; }
      /** reference to underlying H5File object */
      const std::string fileName(){ return _fileName; }
      /** reference to access flag */
      const int accessFlag(){ return _accessFlag; }

    private:
      H5File*     _h5File;
      std::string _fileName;
      int         _accessFlag;

  }; // H5IO class

} // namespace agnos
  

#endif // IO_HANDLER_H
