
#ifndef PHYSICS_CHANNEL_FLOW_H
#define PHYSICS_CHANNEL_FLOW_H

#include "agnosDefines.h"
#include "PhysicsLibmesh.h"

/** channelflow includes */
#include "channel_solver.h"
#include "channel_system.h"


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
      virtual ~PhysicsChannelFlow( );
    
    protected:
      /** set parameter values */
      virtual void _setParameterValues( const T_S& parameterValues ) ;

      /** Redefine _solve( ) routine to call rely on _flowSolver */
      void _solve( );

      /** ChannelSolver object */
      std::shared_ptr<ChannelSolver> _flowSolver ;
      


  };

}

#endif // PHYSICS_CHANNEL_FLOW_H
