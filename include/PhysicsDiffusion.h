

#ifndef PHYSICS_DIFFUSION_H
#define PHYSICS_DIFFUSION_H

#include "agnosDefines.h"
#include "PhysicsLibmesh.h"

// libmesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
        
#define _USE_MATH_DEFINES



namespace AGNOS
{



  /********************************************//**
   * \brief Example PhysicsLibmesh class - higher order diffusion example from
   * paper
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsDiffusion : public PhysicsLibmesh<T_S,T_P>
  {

    public:
      /** Constructor. Pass input file to provide setting to physics class */
      PhysicsDiffusion( const Communicator& comm_in, const GetPot& input );

      /** destructor */
      virtual ~PhysicsDiffusion( );


    protected:
      /** Geometry and boundary data */
      unsigned int    _nElem;
      
      /** set parameter values */
      virtual void _setParameterValues( const T_S& parameterValues ) ;

  };



}
#endif // PHYSICS_DIFFUSION_H
