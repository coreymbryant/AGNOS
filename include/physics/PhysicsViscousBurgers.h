
#ifndef PHYSICS_VISCOUS_BURGERS_H
#define PHYSICS_VISCOUS_BURGERS_H


#include "agnosDefines.h"
#include "PhysicsLibmesh.h"

// libmesh includes
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/petsc_nonlinear_solver.h"
#include "libmesh/petsc_diff_solver.h"
#include "libmesh/newton_solver.h"
#include "libmesh/steady_solver.h"

namespace AGNOS
{

  /********************************************//**
   * \brief Basic 1D Burger's PhysicsModel class
   *
   * This example is given in the book
   * "Spectral Methods for Uncertainty Quantification" by Le Maitre and Knio
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsViscousBurgers : public PhysicsLibmesh<T_S,T_P>
  {

    public:
      /** Constructor. Pass input file to provide setting to physics class */
      PhysicsViscousBurgers( const Communicator& comm_in, const GetPot& input );

      /** Destructor */
      virtual ~PhysicsViscousBurgers( );

      /** Redefine exactQoi for this model */
      T_P exactQoi( )
      {
        T_P resultVector(1);
        resultVector(0) = 10.;
        return resultVector;
      }

    protected:
      /** Geometry and boundary data */
      double          _L;
      int             _nElem;

      Number          _mu;

      /** solver settings */ 
      unsigned int    _nonlinearTolerance;
      unsigned int    _nonlinearSteps;


      /** set parameter values */
      virtual void _setParameterValues( const T_S& parameterValues ) ;

      

  };




}

#endif // PHYSICS_VISCOUS_BURGERS_H
