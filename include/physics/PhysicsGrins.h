
#ifndef PHYSICS_GRINS_FLOW_H
#define PHYSICS_GRINS_FLOW_H

#include <fstream>

#include "PhysicsLibmesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmv_io.h"

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

      /** Reference to qoi value */
      std::shared_ptr<libMesh::DifferentiableQoI> _differentiableQoI ;
      virtual libMesh::Number _getQoIValue( ) 
      { return this->_system->qoi[0]; }

      /** grins input file  */
      GetPot _grinsInput ;

      /** reference for a physicsList */
      GRINS::PhysicsList _physicsList ;

      /** Parameter names and sizes */
      std::map<std::string, std::map<std::string,std::string> > 
        _grinsParameters;

      /** substitute variables with parameter values in input file */
      void _substituteVariable( std::string varName, std::string varValue );

      /** reference to GRINS multiphysics system */
      GRINS::MultiphysicsSystem* _multiphysicsSystem;

      /** Simulation builder */
      std::shared_ptr<GRINS::SimulationBuilder> _simulationBuilder ;

      /** Simulation  */
      std::shared_ptr<GRINS::Simulation> _simulation;

      /** Output file */
      libMesh::ExodusII_IO* _adjointExio;
      libMesh::ExodusII_IO* _primalExio;
      void _initOutput( )
      {
        _adjointExio = new ExodusII_IO(*this->_mesh);
        _primalExio = new ExodusII_IO(*this->_mesh);
        return;
      }

      /** Redefine how primal solution is output */
      void _outputPrimalViz( const T_S& paramVector )
      {
        std::string fileName = "primal.exo" ;
        
        // output timestep info
        _primalExio->write_timestep(
            fileName,
            *this->_equationSystems,
            this->_timeStep, this->_timeStep );

        // output solution vectors
        std::vector<std::string> solNames;
        std::vector<libMesh::Number> solVec;
        this->_equationSystems->build_variable_names(solNames,NULL,NULL);
        this->_equationSystems->build_solution_vector(solVec,NULL);
        _primalExio->write_nodal_data(
            fileName,
            solVec,
            solNames
            );

        // output parameters as global vars
        std::vector<std::string> paramNames;
        for(unsigned int i=0; i<paramVector.size();i++)
          paramNames.push_back("parameter"+std::to_string(i));

        std::vector<Number> soln = paramVector.get_values();

        _primalExio->write_global_data( soln, paramNames );

        return;
      }

      /** Redefine how adjoint solution is output */
      void _outputAdjointViz( const T_S& paramVector )
      {
        std::string fileName = "adjoint.exo";

        NumericVector<Number> &primal_solution = *this->_system->solution;
        NumericVector<Number> &dual_solution 
          = this->_system->get_adjoint_solution(0);
        primal_solution.swap(dual_solution);

        
        // output timestep info
        _adjointExio->write_timestep(
            fileName,
            *this->_equationSystems ,
            this->_timeStep, this->_timeStep );
        
        // output solution vectors
        std::vector<std::string> solNames;
        std::vector<libMesh::Number> solVec;
        this->_equationSystems->build_variable_names(solNames,NULL,NULL);
        this->_equationSystems->build_solution_vector(solVec,NULL);
        _adjointExio->write_nodal_data(
            fileName,
            solVec,
            solNames
            );

        // output parameters as global vars
        std::vector<std::string> paramNames;
        for(unsigned int i=0; i<paramVector.size();i++)
          paramNames.push_back("parameter"+std::to_string(i));

        std::vector<Number> soln = paramVector.get_values();

        _adjointExio->write_global_data( soln, paramNames );

        primal_solution.swap(dual_solution);
      }

  };



  class MyQoI : public libMesh::DifferentiableQoI
  {
    public:
    /* MyQoI( std::string& qoi_name ) : QoIBase(qoi_name){ return; } */ 
    /* bool assemble_on_interior() const { return true; } */
    /* bool assemble_on_sides() const { return false; } */
      libMesh::AutoPtr<libMesh::DifferentiableQoI> clone() 
      { return libMesh::AutoPtr<libMesh::DifferentiableQoI>(new MyQoI( *this ) ); }
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
      const std::vector<Point > &xyz = elem_fe->get_xyz();

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
        Point x = xyz[qp] ;
        Number u = c.interior_value(0, qp);
        // 10/pi*exp(-10*(x-1).^2 - 10*(y-0).^2) 
        Q[0] += JxW[qp] * (
            u * 10./libMesh::pi * 
            std::exp( -10. * std::pow(x(0)-1.,2.) - 10.*std::pow(x(1)-0.,2.) )
            );
        // u
        /* Q[0] += JxW[qp] * u  ; */
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
      const std::vector<Point > &xyz = elem_fe->get_xyz();

      // The number of local degrees of freedom in each variable
      const unsigned int n_u_dofs = c.get_dof_indices(0).size();
      unsigned int n_qpoints = c.get_element_qrule().n_points();

      // Fill the QoI RHS corresponding to this QoI. Since this is the 0th QoI
      // we fill in the [0][i] subderivatives, i corresponding to the variable index.
      // Our system has only one variable, so we only have to fill the [0][0] subderivative
      DenseSubVector<Number> &Q = c.get_qoi_derivatives(0,0);

      // Loop over the qps
      for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        Point x = xyz[qp] ;
        for (unsigned int i=0; i != n_u_dofs; i++)
          Q(i) += JxW[qp] * (
              phi[i][qp] * 10./libMesh::pi * 
              std::exp( -10. * std::pow(x(0)-1.,2.) - 10.*std::pow(x(1)-0.,2.) )
              );
            
          /* Q(i) += JxW[qp] * phi[i][qp] ; */
      }

    }
  };

} // end of namespace AGNOS

#endif // PHYSICS_GRINS_FLOW_H
