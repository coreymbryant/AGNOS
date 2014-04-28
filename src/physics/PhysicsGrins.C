
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

    // build physics list: use simulation builder
    _physicsList = 
      _simulationBuilder->build_physics( this->_grinsInput ) ;
    
    // read in parameters unique to this model
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre GRINS specific input " << std::endl;
    _grinsParameters.clear();
    // need a non const GetPot object to search
    GetPot modInput( this->_input );
    for( GRINS::PhysicsListIter it = _physicsList.begin(); 
        it != _physicsList.end(); 
        it++)
    {
      std::map<std::string, std::string> physicsParams; 

      std::string prefixName = "physics/"+it->first+"/" ;
      modInput.set_prefix(prefixName.data()) ;

      std::vector<std::string> paramNames = modInput.get_variable_names() ;
      for(unsigned int p=0; p<paramNames.size();p++)
      {
        std::string value = modInput( paramNames[p], "") ;
        physicsParams.insert( 
              std::pair<std::string,std::string>( paramNames[p], value )
              );

        std::cout << "Init: "<< paramNames[p] << ": " << value << std::endl;
      }
      _grinsParameters.insert( 
          std::pair<std::string, std::map<std::string,std::string> >(
            it->first, physicsParams )
          );

    } // end iterate through physics list


    // initialize simulation 
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre GRINS::Simulation init " << std::endl;
    _simulation.reset( 
        new GRINS::Simulation( 
          this->_grinsInput, *_simulationBuilder, this->_communicator ) 
        );

    //  Get pointers to members of GrinsSystem 
    this->_equationSystems = (_simulation->get_equation_system()).get();
    this->_multiphysicsSystem = &( 
        this->_equationSystems->template
        get_system<GRINS::MultiphysicsSystem>("GRINS")
        );
    this->_system = this->_multiphysicsSystem ;

    // build mesh: use simulation builder
    this->_mesh = &(this->_equationSystems->get_mesh() ) ;


    //  Build mesh refinement object 
    this->_buildMeshRefinement();

    this->_equationSystems->reinit();
    
    // Set up QoIs 
    this->_qois = new libMesh::QoISet;
    std::vector<unsigned int> qoi_indices;
    qoi_indices.push_back(0);
    this->_qois->add_indices(qoi_indices);
    /* this->_qois->set_weight(0, 1.0); */
    _differentiableQoI.reset( new GRINS::MyQoI ) ;
    /* std::string name = "mine" ; */
    /* GRINS::MyQoI myqoi( name ); */
    /* _differentiableQoI->add_qoi( myqoi ) ; */
    /* _differentiableQoI->init( _grinsInput, *this->_multiphysicsSystem ) ; */
    this->_multiphysicsSystem->qoi.resize(1);
    this->_multiphysicsSystem->attach_qoi( _differentiableQoI.get() );
    std::cout << "qoi: " << this->_multiphysicsSystem->qoi.size() << std::endl;

    // Build estimator object 
    this->_buildErrorEstimator();
    
    // we need to run a solve routine on all procs to make sure residuals are
    // set up correctly when we go back to compute residuals with surrogate
    // evaluations
    // TODO: do we need to do this for this class?
    /* this->_system->solve( ); */

    // print out some simulation info
    _simulation->print_sim_info();

    _initOutput() ;

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
    std::map<std::string,std::map<std::string,std::string> >::iterator it;
    for( it = _grinsParameters.begin(); it != _grinsParameters.end();
        it++)
    {
      std::map<std::string, std::string>::iterator pit
        = _grinsParameters[it->first].begin(); 
      for( ; pit != _grinsParameters[it->first].end(); pit++)
      {
        if (AGNOS_DEBUG)
          std::cout << "pre -- " << pit->first << ": " << pit->second 
            << std::endl;

        // create copies of input strings to maniuplate
        std::string values( pit->second );
        std::size_t foundBegin=0, foundEnd=0; 
        foundBegin = values.find("$(",foundEnd) ;
        foundEnd = values.find(")",foundBegin+1) ;
        // while they still exist in the line replace them
        while( foundBegin < std::string::npos )
        {
          // extract parameter component
          std::string var =
            values.substr(foundBegin,foundEnd-foundBegin+1) ;
          unsigned int varNum = std::stoi( var.substr(2,var.length()-3) );

          // make sure given index is less than total num of params
          agnos_assert( ( varNum < parameterValues.size() ) );

          // replace with parameter values
          values.replace( 
              foundBegin,
              var.length(),
              std::to_string( parameterValues(varNum) ) 
              ) ;

          if (AGNOS_DEBUG)
          {
            std::cout << " found begin: " << foundBegin << std::endl;
            std::cout << " found end: " << foundEnd << std::endl;
            std::cout << " var: " << var<< std::endl;
            std::cout << " varNum: " << varNum << std::endl;
          }
          
          // find parameters in variable string
          foundBegin = values.find("$(",foundEnd) ;
          foundEnd = values.find(")",foundBegin+1) ;

        } 

        /* if (AGNOS_DEBUG) */
        {
          std::cout << " post -- " 
            << pit->first << ": " << pit->second << std::endl;
          std::cout << " values -- " 
            << pit->first << ": " << values << std::endl;
        }

        // set grins input to new values
        this->_grinsInput.set( 
            "Physics/"+it->first+"/"+pit->first, 
            values
            ) ;
      }

    }
    

    // construct new physics list based on new parameter values
    //    we have to construct a new list in order to set the value of
    //    parameters that are private in GRINS::Physics classes
    _physicsList = 
      this->_simulationBuilder->build_physics( this->_grinsInput ) ;

    // init variables for new physics
    GRINS::PhysicsListIter git ;
    for(  git = _physicsList.begin(); git != _physicsList.end(); git++)
      git->second->init_variables( _multiphysicsSystem );
    
    // attach new physics list to system
    _multiphysicsSystem->attach_physics_list( _physicsList );


    // reinit the equation system just to make sure all routines will use new
    // values
    this->_equationSystems->reinit() ;

  }


  // Define template types
  template class
    PhysicsGrins<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;
}


