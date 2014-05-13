
#include "ioHandler.h"
#include "EvaluatorPseudoSpectral.h"

namespace AGNOS{

  /********************************************//**
   * \brief Constructor
   *
   * by default it opens fileName with H5F_ACC_EXCL access flag and H5P_DEFAULT
   * create and access propery lists. 
   *
   * 
   ***********************************************/
  H5IO::H5IO( std::string fileName, int accessFlag )
    :
    _h5File( new H5File(fileName, accessFlag, H5P_DEFAULT, H5P_DEFAULT )),
    _fileName(fileName),
    _accessFlag(accessFlag)
  {
  }

  
  H5IO::~H5IO( )
  {
    delete _h5File;
  }

  void H5IO::close()
  {
    _h5File->close();
    return;
  }

  /** write user block to file for reference purposes */
  void H5IO::writeUserBlock(  )
  { return ; }

  /** read in user block and get a reference to H5File object */
  void H5IO::readUserBlock(  )
  { return ; }

  /** write parameter info: dim, type, bounds; to HDF5 file*/
  void H5IO::writeParameters(
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

    delete dataset;
    delete dataspace;

    return;
  }

  /** read parameter info: dim, type, bounds; from HDF5 file*/
  void H5IO::readParameters(
      CommonFG* common,
      std::vector<std::shared_ptr<AGNOS::Parameter> >& parameters)
  {
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
  void H5IO::writeSurrogate( 
      CommonFG* common,
      std::shared_ptr<AGNOS::SurrogateModelBase<T_S,T_P> >& surrogate ) 
  {
    int rank = 1;
    hsize_t dims[1];

    DataSpace* dataspace;
    DataSet* dataset;

    // write surrogate parameters data
    Group* group = new Group( common->createGroup( "parameters" ) ) ;
    std::vector<std::shared_ptr<AGNOS::Parameter> > parameters
      = surrogate->getParameters( ) ;
    writeParameters( group, parameters ) ;
    delete group;

    // order data
    std::vector<unsigned int>  order = surrogate->getExpansionOrder();
    dims[0] = order.size();
    dataspace = new DataSpace(rank,dims);
    dataset = new DataSet(
        common->createDataSet("order", PredType::NATIVE_INT, *dataspace ) 
        );
    dataset->write( order.data(), PredType::NATIVE_INT );
    delete dataset;
    delete dataspace;


    // compute solutions names
    std::set<std::string> computeSolutions = surrogate->getSolutionNames();
    dims[0] = computeSolutions.size();
    std::vector<const char*> arr_c_str;
    std::set<std::string>::iterator sit = computeSolutions.begin();
    for(;sit!=computeSolutions.end();sit++)
      arr_c_str.push_back(sit->c_str());
    DataType strDataType = StrType(PredType::C_S1, H5T_VARIABLE);
    dataspace = new DataSpace(rank,dims);
    dataset = new DataSet(
        common->createDataSet("computeSolutions", strDataType, *dataspace ) 
        );
    dataset->write( arr_c_str.data(), strDataType );
    delete dataset;
    delete dataspace;

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

    dataspace = new DataSpace(rank,setDims);
    dataset = new DataSet(
        common->createDataSet("indexSet", PredType::NATIVE_INT, *dataspace ) 
        );
    dataset->write( is, PredType::NATIVE_INT );
    delete dataset;
    delete dataspace;
    
    // coefficients
    group = new Group( common->createGroup( "coefficients" ) ) ;
    std::map< std::string, LocalMatrix > coefficients
      = surrogate->getCoefficients();
    std::map< std::string, LocalMatrix >::iterator cit
      = coefficients.begin(); 

    rank = 2;
    setDims[0] = cit->second.m();
    setDims[1] = cit->second.n();
    dataspace = new DataSpace(rank,setDims);

    nComp = setDims[0] * setDims[1] ;
    double *coeff = new double[nComp];
    for(;cit!=coefficients.end();cit++)
    {
      for(unsigned int i=0; i<setDims[0];i++)
        for(unsigned int j=0; j<setDims[1];j++)
          coeff[i*setDims[1]+j] = cit->second(i,j) ;

      dataset = new DataSet(
          group->createDataSet( 
            cit->first, PredType::NATIVE_DOUBLE, *dataspace) 
          );
      dataset->write( coeff, PredType::NATIVE_DOUBLE );
      delete dataset;
    }
    delete dataspace;
    delete group;

    return;
  }

  /** read surrogate properties from HDF5 file */
  void H5IO::readSurrogate( 
      CommonFG* common,
      std::shared_ptr<AGNOS::SurrogateModelBase<T_S,T_P> >& surrogate,
      const Communicator& comm) 
  {
    DataSet dataset ;
    DataSpace dataspace ;
    int rank = 1;
    hsize_t dims[rank];
    

    // read surrogate parameters data
    std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;
    Group group = common->openGroup("parameters");
    readParameters( &group, parameters ) ;

    // order data
    std::vector<unsigned int> order;
    dataset = common->openDataSet("order");
    dataspace = dataset.getSpace();
    rank = dataspace.getSimpleExtentDims( dims, NULL);
    order.resize(dims[0]);
    dataset.read( order.data(), PredType::NATIVE_INT );

    // compute solutions
    std::set<std::string> computeSolutions;
    dataset = common->openDataSet("computeSolutions");
    dataspace = dataset.getSpace();
    rank = dataspace.getSimpleExtentDims( dims, NULL);

    std::vector<const char*> arr_c_str;
    arr_c_str.resize(dims[0]);
    DataType strDataType = StrType(PredType::C_S1, H5T_VARIABLE);

    dataset.read( arr_c_str.data(), strDataType );

    for(unsigned int i=0; i<arr_c_str.size();i++)
      computeSolutions.insert( std::string(arr_c_str[i]) ) ;

    // index_set
    std::vector< std::vector<unsigned int> > indexSet;
    dataset = common->openDataSet("indexSet");
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
    delete is;

    // coefficients
    std::map< std::string, LocalMatrix > coefficients;
    group = common->openGroup("coefficients");
    hsize_t nObj = group.getNumObjs();

    for(unsigned int o=0; o<nObj; o++)
    {
      std::string solName = group.getObjnameByIdx(o) ;

      dataset = group.openDataSet( solName );
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
      delete coeff;
    }

    // construct the surrogate model as an evaluator
    //  - if a different type is required after intialization the necessary data
    //  can be transferred
    AGNOS::EvaluatorPseudoSpectral<T_S,T_P>* evalSurrogate
      = new AGNOS::EvaluatorPseudoSpectral<T_S,T_P>(
      comm,
      parameters, 
      order,
      computeSolutions);
    evalSurrogate->build(indexSet,coefficients);
    
    surrogate.reset(  evalSurrogate ) ;

    return;
  }


  /** write physics data */
  void H5IO::writePhysics(
      CommonFG* common,
      std::shared_ptr<PhysicsModel<T_S,T_P> >&  physics )
  {
    return ;
  }

  /** read physics data */
  void H5IO::readPhysics(
      CommonFG* common,
      std::shared_ptr<PhysicsModel<T_S,T_P> >&  physics,
      const Communicator& physicsComm)
  {
    return;
  }


  /** write a single element's data to HDF5 file */
  void H5IO::writeElement( 
      CommonFG* common,
      std::shared_ptr<AGNOS::Element<T_S,T_P> >& element )
  {

    // write element specific data
    int rank = 1;
    hsize_t dims[1] = {1};
    DataSpace* dataspace = new DataSpace(rank,dims);
    DataSet* dataset = new DataSet(
        common->createDataSet("weight", PredType::NATIVE_DOUBLE, *dataspace)
        );
    double weight[] = { element->weight() } ;
    dataset->write( weight, PredType::NATIVE_DOUBLE );
    delete dataset;
    delete dataspace;
    

    // write element parameters data
    std::vector<std::shared_ptr<AGNOS::Parameter> > parameters
      = element->parameters( ) ;
    writeParameters( common, parameters ) ;

    // write element surrogates 
    Group* group = new Group( common->createGroup( "surrogates" ) ) ;
    std::vector<std::shared_ptr<SurrogateModelBase<T_S,T_P> > > surrogates =
      element->surrogates() ;
    unsigned int nSurrogates = surrogates.size();
    for(unsigned int i=0;i<nSurrogates;i++)
    {
      Group* subGroup 
        = new Group( group->createGroup( "surrogate_"+std::to_string(i) ) ) ;
      writeSurrogate( subGroup, surrogates[i] );
      delete subGroup;
    }
    delete group;

    // write element physics
    group = new Group( common->createGroup( "physics" ) ) ;
    std::shared_ptr<PhysicsModel<T_S,T_P> >  physics = element->physics( ) ;
    writePhysics( group, physics );
    delete group;

    return;
  }

  /** read a single element's data from HDF5 file */
  void H5IO::readElement( 
      CommonFG* common,
      std::shared_ptr<AGNOS::Element<T_S,T_P> >& element,
      const Communicator& comm, const Communicator& physicsComm )
  {
    DataSet dataset ;
    DataSpace dataspace ;
    int rank = 1;
    hsize_t dims[rank];

    // read element specific data
    double weight[1] ;
    dataset = common->openDataSet("weight");
    dataspace = dataset.getSpace();
    rank = dataspace.getSimpleExtentDims( dims, NULL );
    dataset.read( weight, PredType::NATIVE_DOUBLE );
    
    // read element parameters data
    std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;
    readParameters( common, parameters ) ;

    // read element surrogates 
    Group group = common->openGroup("surrogates");
    std::vector<std::shared_ptr<SurrogateModelBase<T_S,T_P> > > surrogates ;

    hsize_t nSurrogates = group.getNumObjs();
    for(unsigned int i=0;i<nSurrogates;i++)
    {
      std::shared_ptr<SurrogateModelBase<T_S,T_P> > surrogate ;
      Group* subGroup 
        = new Group( group.openGroup( "surrogate_"+std::to_string(i) ) ) ;
      readSurrogate( subGroup, surrogate, comm );

      surrogates.push_back( surrogate );
      delete subGroup;
    }

    
    // read element physics
    group = common->openGroup("physics");
    std::shared_ptr<PhysicsModel<T_S,T_P> >  physics;
    readPhysics( &group, physics, physicsComm );

    // construct element
    element.reset(
        new AGNOS::Element<T_S,T_P>( parameters, surrogates, physics, weight[0])
        ) ;

    return;
  }


  /** write simulation data */
  void H5IO::writeSimulation( std::list<AGNOS::Element<T_S,T_P> >& activeElems )
  {
    Group* group = new Group( _h5File->createGroup("Elements") );

    int elemNum = 0;
    std::list<AGNOS::Element<T_S,T_P> >::iterator elit
      = activeElems.begin();
    for(;elit != activeElems.end(); elit++)
    {
      Group* subGroup 
        = new Group( group->createGroup( "Elem_"+std::to_string(elemNum++) ) ) ;
      std::shared_ptr<AGNOS::Element<T_S,T_P> > elemPtr( 
          new AGNOS::Element<T_S,T_P>(*elit) 
          );
      writeElement(subGroup,elemPtr) ;
    }

    return ;
  }

  /** read simulation data */
  void H5IO::readSimulation( 
      std::list<AGNOS::Element<T_S,T_P> >& activeElems,
      const Communicator& comm, const Communicator& physicsComm )
  {
    activeElems.clear();
    Group group = _h5File->openGroup("Elements");

    hsize_t nElems = group.getNumObjs();

    for(unsigned int i=0; i<nElems; i++)
    {
      Group subGroup 
        = group.openGroup( "Elem_"+std::to_string(i) );
      std::shared_ptr<AGNOS::Element<T_S,T_P> > elem ;
      readElement(&subGroup,elem,comm,physicsComm) ;
      activeElems.push_front(*elem) ;
    }

    return ;
  }

} // end namespace
