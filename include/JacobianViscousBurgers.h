#ifndef JACOBIAN_VISCOUS_BURGERS_H
#define JACOBIAN_VISCOUS_BURGERS_H

#include "agnosDefines.h"
#include "PhysicsJacobian.h"

namespace AGNOS{

  /********************************************//**
   * \brief PhysicsJacobian class for viscousBurgers model
   ***********************************************/
  template<class T_S>
    class JacobianViscousBurgers : public PhysicsJacobian<T_S>
  {

    public:
      JacobianViscousBurgers( );

      void jacobian( 
          const NumericVector<Number> &X, 
          SparseMatrix<Number> &J,
          NonlinearImplicitSystem &S);

      void setSystemData( const GetPot& input ) ;

    protected:
      double          m_L;
      double          m_uMinus;
      double          m_uPlus;
  
  };

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    JacobianViscousBurgers<T_S>::JacobianViscousBurgers()
    {}

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    void JacobianViscousBurgers<T_S>::setSystemData(
        const GetPot& input
        )
    { 
      m_L                 = input("physics/L",10.);
      m_uMinus                 = input("physics/uMinus",1.);
      m_uPlus                 = input("physics/uPlus",0.);
      return;
    }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S> 
    void JacobianViscousBurgers<T_S>::jacobian(
        const NumericVector<Number> &X, 
        SparseMatrix<Number> &J,
        NonlinearImplicitSystem &S
        )
  {
    /* std::cout << "test: jacobian begin" << std::endl; */
    // define viscosity based on current parameter values
    double mu = 1. + 0.62 * (this->m_paramValues(0)) 
      + 0.36 * (this->m_paramValues(1))  ;

    EquationSystems &es = S.get_equation_systems();

    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

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
    AutoPtr<QBase> qrule (QBase::build("Gauss",dim, FIFTH) );

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule (qrule.get());

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

    // Define data structures to contain the Jacobian element matrix.
    // Following basic finite element terminology we will denote these
    // "Ke". More detail is in example 3.
    DenseMatrix<Number> Ke;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;

    // Now we will loop over all the active elements in the mesh which
    // are local to this processor.
    // We will compute the element Jacobian contribution.
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

        // Zero the element Jacobian before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).
        Ke.resize (dof_indices.size(),
                   dof_indices.size());

        // Now we will build the element Jacobian.  This involves
        // a double loop to integrate the test funcions (i) against
        // the trial functions (j). Note that the Jacobian depends
        // on the current solution x, which we access using the soln
        // vector.
        //
        for (unsigned int qp=0; qp<qrule->n_points(); qp++)
          {
            Number u = 0;

            for (unsigned int i=0; i<phi.size(); i++)
              u      += phi[i][qp]*X(dof_indices[i]);


            for (unsigned int i=0; i<phi.size(); i++)
              for (unsigned int j=0; j<phi.size(); j++)
                Ke(i,j) += JxW[qp]*(
                    // mu(du_x,v_x)
                      mu * (dphi[j][qp] * dphi[i][qp] )       // mu(du_x,v_x)
                    // -1/2 (du - 2 u du , v_x)
                    - 0.5 * phi[j][qp] * (1. - 2. * u) * dphi[i][qp](0) 
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
            Ke(s,s) += penalty ;
        }

        dof_map.constrain_element_matrix (Ke, dof_indices);

        // Add the element matrix to the system Jacobian.
        J.add_matrix (Ke, dof_indices);
      }

    // That's it.
  }





}
#endif // JACOBIAN_VISCOUS_BURGERS_H


