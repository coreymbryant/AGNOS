#ifndef RESIDUAL_VISCOUS_BURGERS_H
#define RESIDUAL_VISCOUS_BURGERS_H

#include "agnosDefines.h"
#include "PhysicsResidual.h"

namespace AGNOS{

  /********************************************//**
   * \brief PhysicsResidual class for viscousBurgers model
   ***********************************************/
  template<class T_S>
    class ResidualViscousBurgers : public PhysicsResidual<T_S>
  {

    public:
      ResidualViscousBurgers( );

      void residual( 
          const NumericVector<Number> &X, 
          NumericVector<Number> &R,
          NonlinearImplicitSystem &S) ;

      void setSystemData( const GetPot& input ) ;

    protected:
      double m_min;
      double m_max;
      
  
  };

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    ResidualViscousBurgers<T_S>::ResidualViscousBurgers()
    { }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    void ResidualViscousBurgers<T_S>::setSystemData(
        const GetPot& input
        )
    { 
      m_min = input("physics/min",-10.);
      m_max = input("physics/max",10.);
      return;
    }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S> 
    void ResidualViscousBurgers<T_S>::residual(
        const NumericVector<Number> &X, 
        NumericVector<Number> &R,
        NonlinearImplicitSystem &S
        )
    { 
      /* std::cout << "test: residual begin" << std::endl; */

      // define viscosity based on current parameter values
      double mu = 1. + 0.62 * (this->m_paramValues(0)) 
        + 0.36 * (this->m_paramValues(1))  ;

      EquationSystems &es = S.get_equation_systems();

      // Get a constant reference to the mesh object.
      const MeshBase& mesh = es.get_mesh();

      // The dimension that we are running
      const unsigned int dim = mesh.mesh_dimension();
      libmesh_assert_equal_to (dim, 1);

      // Get a reference to the NonlinearImplicitSystem we are solving
      NonlinearImplicitSystem& system =
        es.get_system<NonlinearImplicitSystem>("Burgers");

      // A reference to the \p DofMap object for this system.  The \p DofMap
      // object handles the index translation from node and element numbers
      // to degree of freedom numbers.  We will talk more about the \p DofMap
      // in future examples.
      const DofMap& dof_map = system.get_dof_map();

      // Get a constant reference to the Finite Element type
      // for the first (and only) variable in the system.
      FEType fe_type = dof_map.variable_type(0);

      // Build a Finite Element object of the specified type.  Since the
      // \p FEBase::build() member dynamically creates memory we will
      // store the object as an \p AutoPtr<FEBase>.  This can be thought
      // of as a pointer that will clean up after itself.
      AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));

      // A 5th order Gauss quadrature rule for numerical integration.
      QGauss qrule (dim, FIFTH);

      // Tell the finite element object to use our quadrature rule.
      fe->attach_quadrature_rule (&qrule);

      // Here we define some references to cell-specific data that
      // will be used to assemble the linear system.
      // We begin with the element Jacobian * quadrature weight at each
      // integration point.
      const std::vector<Real>& JxW = fe->get_JxW();

      // The element shape functions evaluated at the quadrature points.
      const std::vector<std::vector<Real> >& phi = fe->get_phi();

      // The element shape function gradients evaluated at the quadrature
      // points.
      const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

      // Define data structures to contain the resdual contributions
      DenseVector<Number> Re;

      // This vector will hold the degree of freedom indices for
      // the element.  These define where in the global system
      // the element degrees of freedom get mapped.
      std::vector<dof_id_type> dof_indices;

      // Now we will loop over all the active elements in the mesh which
      // are local to this processor.
      // We will compute the element residual.
      R.zero();

      MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
      const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

      for ( ; el != end_el; ++el)
        {
          // Store a pointer to the element we are currently
          // working on.  This allows for nicer syntax later.
          const Elem* elem = *el;

          // Get the degree of freedom indices for the
          // current element.  These define where in the global
          // matrix and right-hand-side this element will
          // contribute to.
          dof_map.dof_indices (elem, dof_indices);

          // Compute the element-specific data for the current
          // element.  This involves computing the location of the
          // quadrature points (q_point) and the shape functions
          // (phi, dphi) for the current element.
          fe->reinit (elem);

          // We use the resize member here because
          // the number of degrees of freedom might have changed from
          // the last element.  Note that this will be the case if the
          // element type is different (i.e. the last element was a
          // triangle, now we are on a quadrilateral).
          Re.resize (dof_indices.size());

          // Now we will build the residual. This involves
          // the construction of the matrix K and multiplication of it
          // with the current solution x. We rearrange this into two loops:
          // In the first, we calculate only the contribution of
          // K_ij*x_j which is independent of the row i. In the second loops,
          // we multiply with the row-dependent part and add it to the element
          // residual.

          for (unsigned int qp=0; qp<qrule.n_points(); qp++)
            {
              Number u = 0;
              Gradient grad_u;

              for (unsigned int j=0; j<phi.size(); j++)
                {
                  u      += phi[j][qp]*X(dof_indices[j]);
                  grad_u += dphi[j][qp]*X(dof_indices[j]);
                }


              for (unsigned int i=0; i<phi.size(); i++)
                Re(i) += JxW[qp]*(
                    // mu(u_x,v_x)
                      mu * (grad_u * dphi[i][qp] )    
                    // -1/2 (u (1-u) , v_x)
                    - 0.5 * u * ( 1. - u ) * dphi[i][qp](0) 
                    );
            }

          // At this point the interior element integration has
          // been completed.  However, we have not yet addressed
          // boundary conditions.
          // Define the penalty parameter used to enforce the BC's
          double penalty = 1.e10;

          // Loop over the sides of this element. For a 1D element, the "sides"
          // are defined as the nodes on each edge of the element, i.e. 1D elements
          // have 2 sides.
          for(unsigned int s=0; s<elem->n_sides(); s++)
          {
            // If this element has a NULL neighbor, then it is on the edge of the
            // mesh and we need to enforce a boundary condition using the penalty
            // method.
            if(elem->neighbor(s) == NULL)
            {
              Point p = elem->point(s);

              // get solution at this point
              // FIXME: technically only true for interpolating basis
              Number u = X(dof_indices[s]);

              // if its left boundary
              if (p < Point((m_min+m_min)/2.0) )
              {
                Re(s) += penalty * ( 
                    u  - (0.5 * ( 1 + std::tanh( m_min / 4. / 1.0) ) )
                    );
              }

              // if its right boundary
              if (elem->point(s) > Point((m_min+m_min)/2.0) )
              {
                Re(s) += penalty * (  
                    u - (0.5 * ( 1 + std::tanh( m_max / 4. / 1.0) ) )
                    );
              }
            }
          }

          dof_map.constrain_element_vector (Re, dof_indices);
          R.add_vector (Re, dof_indices);
        }

      // That's it.
    }





}
#endif // RESIDUAL_VISCOUS_BURGERS_H


