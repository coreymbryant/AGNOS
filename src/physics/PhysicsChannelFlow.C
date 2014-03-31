
#include "PhysicsChannelFlow.h"

/** channelflow includes */
/* #include "channel_solver.h" */
/* #include "channel_system.h" */


/** libMesh includes */

namespace AGNOS
{

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsChannelFlow<T_S,T_P>::PhysicsChannelFlow(
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
    // - this is handled by the ChannelSolver 
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre channel_input " << std::endl;
    const std::string channelInputFile 
      = input("channel_input","./flow.in");

    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre ChannelSolver init " << std::endl;
    _flowSolver.reset( new ChannelSolver(channelInputFile) );

    _flowSolver->get_es().print_info();

    //  Get pointers to members of ChannelSystem 
    this->_equationSystems = &(_flowSolver->get_es());
    this->_system = &(
        this->_equationSystems->template get_system<ChannelSystem>("flow")
        ) ;
    this->_mesh = &(this->_system->get_mesh()); 
    //  Build mesh refinement object 
    this->_buildMeshRefinement();
    // Set up QoIs 
    this->_qois = new libMesh::QoISet;
    std::vector<unsigned int> qoi_indices;
    qoi_indices.push_back(0);
    this->_qois->add_indices(qoi_indices);
    /* this->_qois->set_weight(0, 1.0); */
    // Build estimator object 
    this->_buildErrorEstimator();
    


    // TODO how are we going to take care of output 
    /* _flowSolver->output(); */

    if (AGNOS_DEBUG)
      std::cout << "DEBUG: post ChannelSolver init " << std::endl;

  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsChannelFlow<T_S,T_P>::~PhysicsChannelFlow()
  {

    if (AGNOS_DEBUG)
      std::cout << "DEBUG: post PhysicsChannelFlow destructor " << std::endl;
  }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
  void PhysicsChannelFlow<T_S,T_P>::_solve( )
  {
    // Solve the flow
    try
    {
      _flowSolver->solve( this->_communicator.rank() );
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
  void PhysicsChannelFlow<T_S,T_P>::_setParameterValues(
    const T_S& parameterValues )
  {

    // Convert T_S vector to stl vector before calling turbulence model
    // setParameters()
    const std::vector<double> params = parameterValues.get_values();

    dynamic_cast<ChannelSystem*>(this->_system)->get_turbulence_model().setParameters(
        params);
  }


  template class
    PhysicsChannelFlow<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;

}

