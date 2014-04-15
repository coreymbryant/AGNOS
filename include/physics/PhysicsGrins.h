
#ifndef PHYSICS_GRINS_FLOW_H
#define PHYSICS_GRINS_FLOW_H

#include "PhysicsLibmesh.h"

/** grins includes */
#include "grins/simulation_builder.h"
#include "grins/simulation.h"
#include "grins/multiphysics_sys.h"


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
  class PhysicsGrins : public PhysicsLibmesh<T_S,T_P>
  {

    public:
      /** Constructor. Pass input file to provide setting to physics class */
      PhysicsGrins( const Communicator& comm_in, const GetPot& input );

      /** Destructor */
      virtual ~PhysicsGrins( );
    
    protected:
      /** Set primal solution to provided vector */
      /* virtual void _setPrimalSolution( T_P& solutionVector ); */

      /** set parameter values */
      virtual void _setParameterValues( const T_S& parameterValues ) ;

      /** Redefine _solve( ) routine to call rely on _flowSolver */
      void _solve( );


      /** Simulation builder */
      std::shared_ptr<GRINS::SimulationBuilder> _simulationBuilder ;

      /** Simulation  */
      std::shared_ptr<GRINS::Simulation> _simulation;

  };

}

#endif // PHYSICS_GRINS_FLOW_H
