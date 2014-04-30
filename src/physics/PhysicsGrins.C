
#include "PhysicsGrins.h"
#include "libmesh/dirichlet_boundaries.h"
#include "grins/composite_function.h"

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
    _grinsInput(std::string( input("grins_input","")) )
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
    _differentiableQoI.reset( new MyQoI ) ;
    /* std::string name = "mine" ; */
    /* MyQoI myqoi( name ); */
    /* _differentiableQoI->add_qoi( myqoi ) ; */
    /* _differentiableQoI->init( _grinsInput, *this->_multiphysicsSystem ) ; */
    this->_multiphysicsSystem->qoi.resize(1);
    this->_multiphysicsSystem->attach_qoi( _differentiableQoI.get() );
    std::cout << "qoi: " << this->_multiphysicsSystem->qoi.size() << std::endl;

    // Build estimator object 
    this->_buildErrorEstimator();
    
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
   * \brief set Parameter values for GRINS models
   *
   * We assume that the variables are represented in the standard grins_input
   * file as scalars. Use GetPot bracket expressions if they are not. Also
   * parameters should be listed in the AGNOS input file as 
   *
   * variable_name = '$(variable_number)'
   *
   * This will ensure they are properly replaced when setting parameter values
   * for AGNOS.
   *
   * Will also work if one needs to define a vector variable based on one
   * parameter, e.g. 
   *
   * parabolic_profile_coeffs_1 = '0.0 0.0 $(0) 0.0 0.0 $(1)'
   *
   * but user needs to be careful if evaluations such as $(0)/2.0 are needed
   * since GetPot won't evaluate them for you without bracket expression, but it
   * will try to evaluate bracket expressions on initilization of the GetPot
   * object, before parameter values are inserted. . 
   ***********************************************/
  template<class T_S,class T_P>
  void PhysicsGrins<T_S,T_P>::_setParameterValues(
    const T_S& parameterValues )
  {
    
    // counter for the total number of parameters
    unsigned int nParams = 0;

    // loop through physics in physics list
    std::map<std::string,std::map<std::string,std::string> >::iterator it;
    for( it = _grinsParameters.begin(); it != _grinsParameters.end();
        it++)
    {

      // loop through variable definitions
      std::map<std::string, std::string>::iterator pit
        = _grinsParameters[it->first].begin(); 
      for( ; pit != _grinsParameters[it->first].end(); pit++)
      {
        // output info
        /* if (AGNOS_DEBUG) */
          std::cout << "pre -- " << pit->first << ": " << pit->second 
            << std::endl;

        // create copies of input strings to maniuplate
        std::string values( pit->second );

        // initialize positions to 0
        std::size_t foundBegin=0, foundEnd=0; 

        // find begining and end of parameter insertion
        foundBegin = values.find("$(",foundEnd) ;
        foundEnd = values.find(")",foundBegin+1) ;

        // if and while params still exist in the line replace them
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

        // find and replace variable strings in grins_input before reparsing 
        _substituteVariable(std::string(pit->first),values);

        /* if (AGNOS_DEBUG) */
        {
          std::cout << " post -- " 
            << pit->first << ": " << pit->second << std::endl;
          std::cout << " values -- " 
            << pit->first << ": " << values << std::endl;
        }

      } // end loop over variable names

    } // end loop over physics in physics list

    
    // reparse grins input file to set new values
    this->_grinsInput = GetPot( this->_input("physics/grins_input","") );

    
    // construct new physics list based on new parameter values
    //    we have to construct a new list in order to set the value of
    //    parameters that are private in GRINS::Physics classes
    _physicsList = 
      this->_simulationBuilder->build_physics( this->_grinsInput ) ;


    // get a reference to dof_map and dirichlet boundaries
    libMesh::DofMap* dof_map = &(_multiphysicsSystem->get_dof_map());
    libMesh::DirichletBoundaries* dirichlet_boundaries 
      = dof_map->get_dirichlet_boundaries();

    // remove all old dirichlet boundaries
    while( !(dof_map->get_dirichlet_boundaries()->empty()) )
      dof_map->remove_dirichlet_boundary(*((*dirichlet_boundaries)[0]));

    // TODO this wont work for periodic boundaries because DofMap doesn't have a
    // remove_periodic_boundary method
    /* libMesh::PeriodicBoundaries* periodic_boundaries */ 
      /* = dof_map->get_periodic_boundaries(); */
    /* while( !(dof_map->get_periodic_boundaries()->empty()) ) */
    /*   dof_map->remove_periodic_boundary(*((*periodic_boundaries)[0])); */

    // reinit the physics variables
    {
      GRINS::PhysicsListIter git = _physicsList.begin();
      for(  ; git != _physicsList.end(); git++ )
        git->second->init_variables( _multiphysicsSystem );

      git = _physicsList.begin();
      for(  ; git != _physicsList.end(); git++ )
        git->second->set_time_evolving_vars( _multiphysicsSystem );

      git = _physicsList.begin();
      git->second->set_is_steady((_multiphysicsSystem->time_solver)->is_steady());

      for(  ; git != _physicsList.end(); git++ )
        git->second->init_bcs( _multiphysicsSystem );

      GRINS::CompositeFunction<libMesh::Number> ic_function;
      for(  ; git != _physicsList.end(); git++ )
        git->second->init_ics( _multiphysicsSystem, ic_function );

      if (ic_function.n_subfunctions())
          _multiphysicsSystem->project_solution(&ic_function);
    }


    // attach new physics list to system
    _multiphysicsSystem->attach_physics_list( _physicsList );
    _multiphysicsSystem->read_input_options( this->_grinsInput );


  }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
  void PhysicsGrins<T_S,T_P>::_substituteVariable(
      std::string varName, std::string varValue )
  {
    // open grins input file
    std::string inputFile = this->_input("physics/grins_input","");
    std::ifstream in( inputFile  ) ;
    std::ofstream out( inputFile+".tmp" , std::ofstream::trunc ) ;

    // vector to hold all of the file
    std::vector<std::string> lines;


    // read in lines
    std::string line;
    while(std::getline(in,line))
    {
      lines.push_back(line);
    }

    //close file
    in.close();

    for(unsigned int l=0;l<lines.size();l++)
    {

      // find variable in input file
      size_t found = lines[l].find(varName+" = ");
      // if found replace value with new value
      if ( found != std::string::npos )
        lines[l].replace(found+varName.length(),std::string::npos,
            " = '"+varValue+"'") ;

    }

    // reopen file and overwrite
    for(unsigned int l=0;l<lines.size();l++)
      out << lines[l] << "\n" ;
    out.close();

    std::rename((inputFile+".tmp").c_str(),inputFile.c_str());

    return; 
  }


  // Define template types
  template class
    PhysicsGrins<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;
}


