
#include "PhysicsGrins.h"

/** libMesh includes */

namespace AGNOS
{

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsGrins<T_S,T_P>::PhysicsGrins(
        const Communicator& comm_in, 
        const GetPot& input 
        ) 
  :
    PhysicsLibmesh<T_S,T_P>(comm_in,input),
    _grinsInput(std::string( input("grins_input","./grins.in")) ) 
  {
    // define available solution names
    this->_availableSolutions.insert("primal");
    this->_availableSolutions.insert("adjoint");
    this->_availableSolutions.insert("qoi");
    this->_availableSolutions.insert("errorEstimate");
    this->_availableSolutions.insert("errorIndicators");


    // pointer to a new simulationBuilder
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre GRINS::SimulationBuilder init " << std::endl;
    _simulationBuilder.reset( new GRINS::SimulationBuilder );
    GRINS::PhysicsList physics_list = 
      _simulationBuilder->build_physics( this->_grinsInput ) ;

    // initialize simulation 
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre GRINS::Simulation init " << std::endl;
    _simulation.reset( 
        new GRINS::Simulation( 
          this->_grinsInput, *_simulationBuilder, this->_communicator ) 
        );

    // print out some simulation info
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre print Simulation info " << std::endl;
    _simulation->print_sim_info();
    _simulation->get_equation_system()->print_info();
    
    // read in parameters unique to this model
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre GRINS specific input " << std::endl;
    _parameterNames.clear();
    for( GRINS::PhysicsListIter it = physics_list.begin(); 
        it != physics_list.end(); 
        it++)
    {
      std::string varName = it->first+"/parameterNames" ;

      int nParamNames = this->_input.vector_variable_size(varName);
      std::vector<std::string> paramNames; 
      paramNames.reserve(nParamNames);

      for(unsigned int p=0; p<nParamNames; p++)
        paramNames.push_back( this->_input(varName,"",p) ) ;

      _parameterNames.insert( 
          std::pair<std::string, std::vector<std::string> >(
            it->first, paramNames )
          );

    }


    //  Get pointers to members of GrinsSystem 
    this->_equationSystems = (_simulation->get_equation_system()).get();
    this->_system = &(
        this->_equationSystems->template
        get_system<GRINS::MultiphysicsSystem>("GRINS")
        ) ;
    this->_multiphysicsSystem = ( 
      static_cast<GRINS::MultiphysicsSystem*>(this->_system)
        );
    this->_mesh = &(this->_system->get_mesh()); 
    //  Build mesh refinement object 
    this->_buildMeshRefinement();
    // Set up QoIs 
    std::shared_ptr<GRINS::CompositeQoI> qoi( new GRINS::CompositeQoI ) ;
    std::string name = "mine" ;
    GRINS::MyQoI myqoi( name );
    qoi->add_qoi( myqoi ) ;
    qoi->init( _grinsInput, *this->_multiphysicsSystem ) ;
    this->_multiphysicsSystem->attach_qoi( qoi.get() );

    this->_qois = new libMesh::QoISet;
    std::vector<unsigned int> qoi_indices;
    qoi_indices.push_back(0);
    this->_qois->add_indices(qoi_indices);
    /* this->_qois->set_weight(0, 1.0); */
    // Build estimator object 
    this->_buildErrorEstimator();
    
    /* // we need to run a solve routine on all procs to make sure residuals are */
    /* // set up correctly when we go back to compute residuals with surrogate */
    /* // evaluations */
    this->_system->solve( );

  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsGrins<T_S,T_P>::~PhysicsGrins()
  {

    if (AGNOS_DEBUG)
      std::cout << "DEBUG: post PhysicsGrins destructor " << std::endl;
  }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
  void PhysicsGrins<T_S,T_P>::_solve( )
  {
    // Solve the flow
    try
    {
      this->_system->solve( );
    }
    catch(int err)
    {
      // only throw when we successfully run requested number of iters and fail
      // to converge!
      std::cout 
        << "*************** Flow solver failed to converge ****************" 
        << std::endl;
      exit(1);
    }

  }


  /********************************************//**
   * All parameters bust be set. If only a subset of parameters are being
   * treated as uncertain then set the deterministic parameters to CONSTANT.
   * This will cause them to be treated as deterministic even though they will
   * still be seen by AGNOS. 
   *
   * NOTE: This may cause issues in the future if we want try to do anisotropic
   * refinement. 
   *
   ***********************************************/
  template<class T_S,class T_P>
  void PhysicsGrins<T_S,T_P>::_setParameterValues(
    const T_S& parameterValues )
  {

    // change input file values 
    //  - only set parameters that the user sets in in grins_input
    //  - this way any type of physics can be handled from this class
    unsigned int nParams = 0;
    std::map<std::string,std::vector<std::string> >::iterator it;
    for( it = _parameterNames.begin(); 
        it != _parameterNames.end(); 
        it++)
      for(unsigned int p=0; p<it->second.size(); p++,nParams++)
        this->_grinsInput.set( 
            "Physics/"+it->first+"/"+it->second[p], 
            parameterValues(p)
            ) ;
    
    // make sure provided parameterValue size matches physics variables size
    // i.e. that number of AGNOS::Paramaters is the same as
    // parameterNames.size()
    agnos_assert( ( parameterValues.size() == nParams ) );

    // construct new physics list based on these parameter values
    //    - this is necessary to set the value of parameters that are private in
    //    GRINS::Physics classes
    GRINS::PhysicsList physics_list = 
      this->_simulationBuilder->build_physics( this->_grinsInput ) ;
    this->_multiphysicsSystem->attach_physics_list( physics_list );
    this->_multiphysicsSystem->read_input_options( this->_grinsInput );
    // have to call init again because newly attached physics_list creates a new
    // uninitialied system ??
    this->_equationSystems->init() ;

  }


  // Define template types
  template class
    PhysicsGrins<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;
}


