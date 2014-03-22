
#ifndef PHYSICS_CHANNEL_FLOW_H
#define PHYSICS_CHANNEL_FLOW_H

#include "agnosDefines.h"
#include "PhysicsLibmesh.h"

/** channelflow includes */
#include "channel_solver.h"


/** libMesh includes */

namespace AGNOS
{


  /********************************************//**
   * \brief PhysicsModel class for running turbulence simulations using
   * channelflow code.
   *
   * 
   ***********************************************/

  template<class T_S, class T_P>
  class PhysicsChannelFlow : public PhysicsLibmesh<T_S,T_P>
  {

    public:
      /** Constructor. Pass input file to provide setting to physics class */
      PhysicsChannelFlow( const Communicator& comm_in, const GetPot& input );

      /** Destructor */
      /* virtual ~PhysicsChannelFlow( ); */
    
    protected:
      /** set parameter values */
      virtual void _setParameterValues( const T_S& parameterValues ) ;

      /** ChannelSolver object */
      std::shared_ptr<ChannelSolver> flowSolver ;
      


  };

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
    // read in parameters unique to this model
    // -> this is handled by the ChannelSolver 
    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre channel_input " << std::endl;
    const std::string channelInputFile 
      = input("channel_input","./flow.in");

    if (AGNOS_DEBUG)
      std::cout << "DEBUG: pre ChannelSolver init " << std::endl;
    flowSolver.reset( new ChannelSolver(channelInputFile) );

    flowSolver->get_es().print_info();


    /* // Get parameter vector from input file */
    /* std::vector<double>* params=NULL; */
    /* std::vector<double> tmp; */
    /* if (flowSolver->myControl.have_variable("ModelParams")) */
    /* { */
    /*   tmp = flowSolver->myControl.vector_value<double>("ModelParams"); */
    /*   params = &tmp; */
    /* } */
    
    /* // Solve the flow */
    /* try */
    /* { */
    /*   flowSolver->solve(0, params); */
    /* } */
    /* catch(int err) */
    /* { */
    /*   // only throw when we successfully run requested number of iters and fail */
    /*   // to converge! */
    /*   std::cout */ 
    /*     << "*************** Flow solver failed to converge ****************" */ 
    /*     << std::endl; */
    /*   exit(1); */
    /* } */


    /* flowSolver->output(); */

    if (AGNOS_DEBUG)
      std::cout << "DEBUG: post ChannelSolver init " << std::endl;

  }

/********************************************//**
 * \brief 
 ***********************************************/
  /* template<class T_S, class T_P> */
  /* PhysicsChannelFlow<T_S,T_P>::~PhysicsChannelFlow() */
  /* { */
  /*   /1* delete flowSolver ; *1/ */

  /*   if (AGNOS_DEBUG) */
  /*     std::cout << "DEBUG: post PhysicsChannelFlow destructor " << std::endl; */

  /* } */

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
  void PhysicsChannelFlow<T_S,T_P>::_setParameterValues(
    const T_S& parameterValues )
  {
    // TODO
    /* static_cast<BurgersSystem*>(this->_system)->_mu */ 
    /*   = 1.0 + 0.62 * parameterValues(0) + 0.36 * parameterValues(1) ; */
  }


}

#endif // PHYSICS_CHANNEL_FLOW_H
