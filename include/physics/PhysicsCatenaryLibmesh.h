

#ifndef PHYSICS_CATENARY_LIBMESH_H
#define PHYSICS_CATENARY_LIBMESH_H


#include "agnosDefines.h"
#include "PhysicsLibmesh.h"

// libmesh includes
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/steady_solver.h"
#include "libmesh/newton_solver.h"
        



namespace AGNOS
{

  /********************************************//**
   * \brief Example PhysicsLibmesh class - catenary chain (solution computed
   * using libmesh)
   *
   * A simple 1D example useful for testing purposes
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsCatenaryLibmesh : public PhysicsLibmesh<T_S,T_P>
  {

    public:
      /** Constructor. Pass input file to provide setting to physics class */
      PhysicsCatenaryLibmesh( const Communicator& comm_in, const GetPot& input );

      /** destructor */
      virtual ~PhysicsCatenaryLibmesh( );

      /** Redefine exactQoi for this model */
      T_P exactQoi( ) ;

    protected:
      /** Geometry and boundary data */
      double          _min;
      double          _max;
      unsigned int    _nElem;

      double          _forcing;

      /** set parameter values */
      virtual void _setParameterValues( const T_S& parameterValues ) ;

  };

}

#endif // PHYSICS_CATENARY_LIBMESH_H
