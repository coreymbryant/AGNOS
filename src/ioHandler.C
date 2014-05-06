
#include "ioHandler.h"

namespace AGNOS{

  /** write user block to file for reference purposes */
  H5File* writeUserBlock( std::string fileName )
  {
    // create a new file using default properties
    H5File* h5_file 
      = new H5File(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    return h5_file;
  }

  /** read in user block and get a reference to H5File object */
  // H5F_RDONLY
  H5File* readUserBlock( std::string fileName )
  {
    // create a new file using default properties
    H5File* h5_file 
      = new H5File(fileName, H5F_ACC_RDONLY, H5P_DEFAULT, H5P_DEFAULT );

    return h5_file;
  }

  /** write parameter info: dim, type, bounds; to HDF5 file*/
  void writeParameters(
      CommonFG* common,
      std::vector<std::shared_ptr<AGNOS::Parameter> > parameters)
  {

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
    int rank = 1;
    hsize_t dims[1] = { parameters.size() } ;
    DataSpace* dataspace = new DataSpace( rank, dims );

    // create composite datatype
    CompType mtype( sizeof(param_type) );
    mtype.insertMember("type",HOFFSET(param_type, type), PredType::NATIVE_INT);
    mtype.insertMember("min",HOFFSET(param_type, min), PredType::NATIVE_DOUBLE);
    mtype.insertMember("max",HOFFSET(param_type, max), PredType::NATIVE_DOUBLE);

    DataSet* dataset = new DataSet( 
        common->createDataSet("parameters",
          mtype, *dataspace) ) ;
    dataset->write( p, mtype );

    delete dataspace;
    delete dataset;

    return;
  }

  /** read parameter info: dim, type, bounds; from HDF5 file*/
  void readParameters(
      CommonFG* common,
      std::vector<std::shared_ptr<AGNOS::Parameter> >& parameters)
  {
    // TODO
    /* Group group = common->openGroup("parameters"); */
    //surrogate model 

    DataSet dataset = common->openDataSet("parameters");
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    hsize_t dims[rank];
    rank = dataspace.getSimpleExtentDims( dims, NULL);
    
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
    param_type p[dims[0]];
    dataset.read(p,mtype);
    parameters.clear();
    for(unsigned int i=0;i<dims[0];i++)
      parameters.push_back(
          std::shared_ptr<AGNOS::Parameter>(
            new AGNOS::Parameter( p[i].type, p[i].min, p[i].max ) ) 
          );

    return;
  }

  /** write surrogate properties to HDF5 file */
  void writeSurrogate( 
      CommonFG* common,
      std::shared_ptr<AGNOS::SurrogateModel<T_S,T_P> >& surrogate )
  {
    /* H5File h5_file(fileName, H5F_ACC_RDWR, H5P_DEFAULT, H5P_DEFAULT ); */

    Group* group = new Group( common->createGroup( "surrogate" ) );

    int rank = 1;
    hsize_t dims[1];

    DataSpace* dataspace;
    DataSet* dataset;

    // order data
    std::vector<unsigned int>  order = surrogate->getExpansionOrder();
    dims[0] = order.size();
    delete dataspace;
    dataspace = new DataSpace(rank,dims);
    delete dataset;
    dataset = new DataSet(
        group->createDataSet("order", PredType::NATIVE_INT, *dataspace ) 
        );
    dataset->write( order.data(), PredType::NATIVE_INT );


    // compute solutions names
    std::set<std::string> computeSolutions = surrogate->getSolutionNames();
    dims[0] = computeSolutions.size();
    std::vector<const char*> arr_c_str;
    std::set<std::string>::iterator sit = computeSolutions.begin();
    for(;sit!=computeSolutions.end();sit++)
      arr_c_str.push_back(sit->c_str());
    DataType strDataType = StrType(PredType::C_S1, H5T_VARIABLE);
    delete dataspace;
    dataspace = new DataSpace(rank,dims);
    delete dataset;
    dataset = new DataSet(
        group->createDataSet("computeSolutions", strDataType, *dataspace ) 
        );
    dataset->write( arr_c_str.data(), strDataType );

    // index_set
    std::vector< std::vector<unsigned int> > indexSet 
      = surrogate->indexSet();
    rank = 2;
    hsize_t setDims[2];
    setDims[0] = indexSet.size(); 
    setDims[1] = indexSet[0].size();
    int nComp = setDims[0] * setDims[1] ;
    int *is = new int[nComp];
    for(unsigned int i=0; i<setDims[0];i++)
      for(unsigned int j=0; j<setDims[1];j++)
        is[i*setDims[1]+j] = indexSet[i][j] ;

    delete dataspace;
    dataspace = new DataSpace(rank,setDims);
    dataset = new DataSet(
        group->createDataSet("indexSet", PredType::NATIVE_INT, *dataspace ) 
        );
    dataset->write( is, PredType::NATIVE_INT );
    
    // coefficients
    Group* subGroup = new Group( group->createGroup( "coefficients" ) ) ;
    std::map< std::string, LocalMatrix > coefficients
      = surrogate->getCoefficients();
    std::map< std::string, LocalMatrix >::iterator cit
      = coefficients.begin(); 

    rank = 2;
    setDims[0] = cit->second.m();
    setDims[1] = cit->second.n();
    delete dataspace;
    dataspace = new DataSpace(rank,setDims);

    nComp = setDims[0] * setDims[1] ;
    double *coeff = new double[nComp];
    for(;cit!=coefficients.end();cit++)
    {
      for(unsigned int i=0; i<setDims[0];i++)
        for(unsigned int j=0; j<setDims[1];j++)
          coeff[i*setDims[1]+j] = cit->second(i,j) ;

      delete dataset;
      dataset = new DataSet(
          subGroup->createDataSet( 
            cit->first, PredType::NATIVE_DOUBLE, *dataspace) 
          );
      dataset->write( coeff, PredType::NATIVE_DOUBLE );
    }

    delete subGroup;


    return;
  }

  /** read surrogate properties from HDF5 file */
  void readSurrogate( 
      CommonFG* common,
      std::vector<unsigned int>& order,
      std::set<std::string>& computeSolutions,
      std::vector< std::vector<unsigned int> >& indexSet,
      std::map< std::string, LocalMatrix >& coefficients
      )
  {

    Group group = common->openGroup( "surrogate" );
    DataSet dataset ;
    DataSpace dataspace ;
    int rank = 1;
    hsize_t dims[rank];

    // order data
    order.clear();
    dataset = group.openDataSet("order");
    dataspace = dataset.getSpace();
    rank = dataspace.getSimpleExtentDims( dims, NULL);
    order.resize(dims[0]);
    dataset.read( order.data(), PredType::NATIVE_INT );

    // compute solutions
    computeSolutions.clear();
    dataset = group.openDataSet("computeSolutions");
    dataspace = dataset.getSpace();
    rank = dataspace.getSimpleExtentDims( dims, NULL);

    std::vector<const char*> arr_c_str;
    arr_c_str.resize(dims[0]);
    DataType strDataType = StrType(PredType::C_S1, H5T_VARIABLE);

    dataset.read( arr_c_str.data(), strDataType );

    for(unsigned int i=0; i<arr_c_str.size();i++)
      computeSolutions.insert( std::string(arr_c_str[i]) ) ;

    // index_set
    indexSet.clear();
    dataset = group.openDataSet("indexSet");
    dataspace = dataset.getSpace();
    hsize_t setDims[2];
    rank = dataspace.getSimpleExtentDims( setDims, NULL);

    int nComp = setDims[0] * setDims[1] ;
    int *is = new int[nComp];
    dataset.read( is, PredType::NATIVE_INT );

    indexSet.resize(setDims[0]);
    for(unsigned int i=0; i<setDims[0];i++)
      for(unsigned int j=0; j<setDims[1];j++)
        indexSet[i].push_back(is[i*setDims[1]+j]) ;

    // coefficients
    coefficients.clear();
    Group subGroup = group.openGroup("coefficients");
    hsize_t nObj = subGroup.getNumObjs();

    for(unsigned int o=0; o<nObj; o++)
    {
      std::string solName = subGroup.getObjnameByIdx(o) ;

      dataset = subGroup.openDataSet( solName );
      dataspace = dataset.getSpace();
      rank = dataspace.getSimpleExtentDims( setDims, NULL);
      nComp = setDims[0] * setDims[1] ;
      double *coeff = new double[nComp];
      dataset.read( coeff, PredType::NATIVE_DOUBLE );

      LocalMatrix coeffMat(setDims[0],setDims[1]);
      for(unsigned int i=0;i<setDims[0];i++)
        for(unsigned int j=0;j<setDims[1];j++)
          coeffMat(i,j) = coeff[i*setDims[1]+j] ;

      coefficients.insert( 
          std::pair<std::string,LocalMatrix>(solName,coeffMat)
          ) ;
    }


    return;
  }


} // end namespace
