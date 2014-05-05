
#include "ioHandler.h"

namespace AGNOS{

  /** write user block to file for reference purposes */
  void writeUserBlock( std::string fileName )
  {
    // create a new file using default properties
    H5File h5_file(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    // close file
    h5_file.close();

    return ;
  }

  /** read in user block and get a reference to H5File object */
  // H5F_RDONLY
  void readUserBlock( std::string fileName )
  {
    // create a new file using default properties
    H5File h5_file(fileName, H5F_ACC_RDONLY, H5P_DEFAULT, H5P_DEFAULT );
    // close file
    h5_file.close();

    return ;
  }

  /** write parameter info: dim, type, bounds; to HDF5 file*/
  void writeParameters(
      std::string fileName,
      std::vector<std::shared_ptr<AGNOS::Parameter> > parameters)
  {
    // open file for read and write
    H5File h5_file(fileName, H5F_ACC_RDWR, H5P_DEFAULT, H5P_DEFAULT );

    // create parameters group
    Group* group = new Group( h5_file.createGroup( "/Parameters") ) ;

    // write dimension data for reference
    int rank = 1;
    hsize_t dims[1] = {1};
    DataSpace* dataspace = new DataSpace(rank,dims);
    DataSet* dataset = new DataSet( 
        h5_file.createDataSet("/Parameters/dimension",
          PredType::NATIVE_INT, *dataspace) ) ;
    int dimension[1] = {(int)parameters.size()};
    dataset->write( dimension, PredType::NATIVE_INT );

    // create datatype for parameter info
    typedef struct param_type {
      int type;
      double min;
      double max;
     } param_type ;

    // fill a vector with parameter data
    param_type p[parameters.size()];
    for(unsigned int i=0;i<parameters.size();i++)
    {
      p[i].type = parameters[i]->type();
      p[i].min = parameters[i]->min();
      p[i].max = parameters[i]->max();
    }

    // create data space 
    rank = 1;
    dims[0] = parameters.size();
    delete dataspace;
    dataspace = new DataSpace( rank, dims );

    // create composite datatype
    CompType mtype( sizeof(param_type) );
    mtype.insertMember("type",HOFFSET(param_type, type), PredType::NATIVE_INT);
    mtype.insertMember("min",HOFFSET(param_type, min), PredType::NATIVE_DOUBLE);
    mtype.insertMember("max",HOFFSET(param_type, max), PredType::NATIVE_DOUBLE);

    delete dataset;
    dataset = new DataSet( 
        h5_file.createDataSet("/Parameters/parameters",
          mtype, *dataspace) ) ;
    dataset->write( p, mtype );

    delete dataspace;
    delete dataset;
    delete group;

    h5_file.close();

    return;
  }

  /** read parameter info: dim, type, bounds; from HDF5 file*/
  void readParameters(
      std::string fileName,
      std::vector<std::shared_ptr<AGNOS::Parameter> >& parameters)
  {
    H5File h5_file(fileName, H5F_ACC_RDONLY, H5P_DEFAULT, H5P_DEFAULT );
    Group group = h5_file.openGroup("/Parameters");

    // read in parameter space dimension
    int dim[1];
    DataSet dataset = group.openDataSet("dimension");
    dataset.read(dim, PredType::NATIVE_INT );

    dataset = group.openDataSet("parameters");
    // create datatype for parameter info
    typedef struct param_type {
      int type;
      double min;
      double max;
     } param_type ;

    // create composite datatype
    CompType mtype( sizeof(param_type) );
    mtype.insertMember("type",HOFFSET(param_type, type), PredType::NATIVE_INT);
    mtype.insertMember("min",HOFFSET(param_type, min), PredType::NATIVE_DOUBLE);
    mtype.insertMember("max",HOFFSET(param_type, max), PredType::NATIVE_DOUBLE);

    // fill a vector with parameter data
    param_type p[dim[0]];
    dataset.read(p,mtype);
    parameters.clear();
    for(unsigned int i=0;i<dim[0];i++)
      parameters.push_back(
          std::shared_ptr<AGNOS::Parameter>(
            new AGNOS::Parameter( p[i].type, p[i].min, p[i].max ) ) 
          );

    // close file
    h5_file.close();
    return;
  }

  /** write surrogate properties to HDF5 file */
  void writeSurrogate( 
      H5File h5_file,
      std::shared_ptr<SurrogateModel<T_S,T_P> > surrogate )
  {
    return;
  }

  /** read surrogate properties from HDF5 file */
  // H5F_RDONLY
  void readSurrogate( 
      H5File h5_file,
      std::shared_ptr<SurrogateModel<T_S,T_P> > surrogate )
  {
    return;
  }


} // end namespace
