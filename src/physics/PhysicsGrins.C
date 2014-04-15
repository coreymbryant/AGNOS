
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
    PhysicsLibmesh<T_S,T_P>(comm_in,input)
  {
    // define available solution names
    this->_availableSolutions.insert("primal");
    this->_availableSolutions.insert("adjoint");
    this->_availableSolutions.insert("qoi");
    this->_availableSolutions.insert("errorEstimate");
    this->_availableSolutions.insert("errorIndicators");

    // read in parameters unique to this model
    /* if (AGNOS_DEBUG) */
      std::cout << "DEBUG: pre grins input " << std::endl;
    const std::string grinsInputFilename = input("grins_input","./grins.in");
    GetPot grinsInputFile( grinsInputFilename ) ;

    /* if (AGNOS_DEBUG) */
      std::cout << "DEBUG: pre GRINS::SimulationBuilder init " << std::endl;
    _simulationBuilder.reset( new GRINS::SimulationBuilder );

    /* if (AGNOS_DEBUG) */
      std::cout << "DEBUG: pre GRINS::Simulation init " << std::endl;
    _simulation.reset( 
        new GRINS::Simulation( 
          grinsInputFile, *_simulationBuilder, this->_communicator ) 
        );

    /* /1* if (AGNOS_DEBUG) *1/ */
    /*   std::cout << "DEBUG: pre print Simulation info " << std::endl; */
    /* _simulation->print_sim_info(); */
    /* _simulation->get_equation_system()->print_info(); */

    /* //  Get pointers to members of ChannelSystem */ 
    /* this->_equationSystems = (_simulation->get_equation_system()).get(); */
    /* //TODO */ 
    /* this->_system = &( */
    /*     this->_equationSystems->template */
    /*     get_system<GRINS::MultiphysicsSystem>("grins") */
    /*     ) ; */
    /* this->_mesh = &(this->_system->get_mesh()); */ 
    /* //  Build mesh refinement object */ 
    /* this->_buildMeshRefinement(); */
    /* // Set up QoIs */ 
    /* this->_qois = new libMesh::QoISet; */
    /* std::vector<unsigned int> qoi_indices; */
    /* qoi_indices.push_back(0); */
    /* this->_qois->add_indices(qoi_indices); */
    /* /1* this->_qois->set_weight(0, 1.0); *1/ */
    /* // Build estimator object */ 
    /* this->_buildErrorEstimator(); */
    
    /* // we need to run a solve routine on all procs to make sure residuals are */
    /* // set up correctly when we go back to compute residuals with surrogate */
    /* // evaluations */
    /* this->_system->solve( ); */

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
      //TODO

  }


  template class
    PhysicsGrins<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;

}

