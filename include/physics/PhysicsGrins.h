
#ifndef PHYSICS_GRINS_FLOW_H
#define PHYSICS_GRINS_FLOW_H

#include "PhysicsLibmesh.h"
#include "libmesh/quadrature.h"

/** grins includes */
#include "grins/simulation_builder.h"
#include "grins/simulation.h"
#include "grins/multiphysics_sys.h"
#include "grins/inc_navier_stokes_base.h"


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

      /** grins input file  */
      GetPot _grinsInput ;

      /** reference for a physicsList */
      GRINS::PhysicsList _physicsList ;

      /** Parameter names and sizes */
      std::map<std::string, std::map<std::string,std::string> > 
        _grinsParameters;

      /** reference to GRINS multiphysics system */
      GRINS::MultiphysicsSystem* _multiphysicsSystem;

      /** Simulation builder */
      std::shared_ptr<GRINS::SimulationBuilder> _simulationBuilder ;

      /** Simulation  */
      std::shared_ptr<GRINS::Simulation> _simulation;

  };

}

namespace GRINS 
{
  class MyQoI : public QoIBase
  {
    public:
    MyQoI( std::string& qoi_name ) : QoIBase(qoi_name){ return; } 
    bool assemble_on_interior() const { return true; }
    bool assemble_on_sides() const { return false; }
    QoIBase* clone() const { return new MyQoI( *this ); }
    //! Compute the qoi value for element interiors.
    virtual void element_qoi( libMesh::DiffContext& context,
                              const libMesh::QoISet& qoi_indices )
    {
      FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

      // First we get some references to cell-specific data that
      // will be used to assemble the linear system.
      FEBase* elem_fe = NULL;
      c.get_element_fe( 0, elem_fe );

      // Element Jacobian * quadrature weights for interior integration
      const std::vector<Real> &JxW = elem_fe->get_JxW();

      // The element quadrature points
      const std::vector<Point > &q_point = elem_fe->get_xyz();

      // The number of local degrees of freedom in each variable
      const unsigned int n_u_dofs = c.get_dof_indices(0).size();
      unsigned int n_qpoints = c.get_element_qrule().n_points();

      // Fill the QoI RHS corresponding to this QoI. Since this is the 0th QoI
      // we fill in the [0][i] subderivatives, i corresponding to the variable index.
      // Our system has only one variable, so we only have to fill the [0][0] subderivative
      std::vector<Number> &Q = c.get_qois();

      // Loop over the qps
      for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        Number u = c.interior_value(0, qp);
        Q[0] += JxW[qp] * u  ;
      } // end of the quadrature point qp-loop
    }

    //! Compute the qoi derivative with respect to the solution on element interiors.
    virtual void element_qoi_derivative( libMesh::DiffContext &context,
                                         const libMesh::QoISet &qoi_indices )
    {
      FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

      // First we get some references to cell-specific data that
      // will be used to assemble the linear system.
      FEBase* elem_fe = NULL;
      c.get_element_fe( 0, elem_fe );

      // Element Jacobian * quadrature weights for interior integration
      const std::vector<Real> &JxW = elem_fe->get_JxW();

      // The basis functions for the element
      const std::vector<std::vector<Real> > &phi = elem_fe->get_phi();

      // The element quadrature points
      const std::vector<Point > &q_point = elem_fe->get_xyz();

      // The number of local degrees of freedom in each variable
      const unsigned int n_u_dofs = c.get_dof_indices(0).size();
      unsigned int n_qpoints = c.get_element_qrule().n_points();

      // Fill the QoI RHS corresponding to this QoI. Since this is the 0th QoI
      // we fill in the [0][i] subderivatives, i corresponding to the variable index.
      // Our system has only one variable, so we only have to fill the [0][0] subderivative
      DenseSubVector<Number> &Q = c.get_qoi_derivatives(0,0);

      // Loop over the qps
      for (unsigned int qp=0; qp != n_qpoints; qp++)
        for (unsigned int i=0; i != n_u_dofs; i++)
          Q(i) += JxW[qp] * phi[i][qp] ;

    }
  };

}

#endif // PHYSICS_GRINS_FLOW_H
