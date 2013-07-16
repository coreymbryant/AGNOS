#ifndef ASSEMBLY_VISCOUS_BURGERS_H
#define ASSEMBLY_VISCOUS_BURGERS_H


#include "agnosDefines.h"

#include "PhysicsAssembly.h"

namespace AGNOS
{


  /********************************************//**
   * \brief PhysicsAssembly class for ViscousBurgers problem.
   ***********************************************/
  template<class T_S>
  class AssemblyViscousBurgers : public PhysicsAssembly
  {
    public:
      AssemblyViscousBurgers(
        libMesh::EquationSystems& es, 
        const std::string&        systemName
          );
      ~AssemblyViscousBurgers( );

      void assemble( );

      void setSystemData(  ) ;

    private:

  };

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    AssemblyViscousBurgers<T_S>::AssemblyViscousBurgers( 
        libMesh::EquationSystems& es, 
        const std::string&        systemName
        )
    : PhysicsAssembly<T_S>(es,systemName)
    { 
    }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    void AssemblyViscousBurgers<T_S>::setSystemData()
    { 
      return;
    }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    AssemblyViscousBurgers<T_S>::~AssemblyViscousBurgers( ){ }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    void AssemblyViscousBurgers<T_S>::assemble( )
   
    {
      // TODO

      // reference to mesh
      const libMesh::MeshBase& mesh = this->m_es.get_mesh();
      const unsigned int dim = mesh.mesh_dimension();
      NonlinearImplicitSystem& system 
        = this->m_es.template get_system<NonlinearImplicitSystem>("1D") ;

      // dofs map
      const libMesh::DofMap& dof_map = system.get_dof_map();
      libMesh::FEType fe_type = dof_map.variable_type(0);

      libMesh::AutoPtr<FEBase> fe(FEBase::build(dim, fe_type));

      // quadrature
      libMesh::QGauss qrule(dim,FOURTH);
      fe->attach_quadrature_rule(&qrule);

      // references to jacobians
      const std::vector<libMesh::Real>& JxW = fe->get_JxW();
      // shape functions at quadrature points
      const std::vector<std::vector<libMesh::Real> >& phi = fe->get_phi();
      const std::vector<std::vector<libMesh::RealGradient> >& dphi = fe->get_dphi();

      // declare matrices
      libMesh::DenseMatrix<libMesh::Number> Ke;
      libMesh::DenseVector<libMesh::Number> Fe;
      std::vector<libMesh::dof_id_type> dof_indices;
      
      // loop over all elements
      libMesh::MeshBase::const_element_iterator el     
        = mesh.active_local_elements_begin();
      const libMesh::MeshBase::const_element_iterator el_end 
        = mesh.active_local_elements_end();

      // Note that ++el is preferred to el++ when using loops with iterators
      for( ; el != el_end; ++el)
      {

        // It is convenient to store a pointer to the current element
        const libMesh::Elem* elem = *el;

        // Get the degree of freedom indices for the current element. 
        // These define where in the global matrix and right-hand-side this 
        // element will contribute to.
        dof_map.dof_indices(elem, dof_indices);

        // Compute the element-specific data for the current element. This 
        // involves computing the location of the quadrature points (q_point) 
        // and the shape functions (phi, dphi) for the current element.
        fe->reinit(elem);

        // Store the number of local degrees of freedom contained in this element
        const int n_dofs = dof_indices.size();

        // We resize and zero out Ke and Fe (resize() also clears the matrix and
        // vector). In this example, all elements in the mesh are EDGE3's, so 
        // Ke will always be 3x3, and Fe will always be 3x1. If the mesh contained
        // different element types, then the size of Ke and Fe would change.
        Ke.resize(n_dofs, n_dofs);
        Fe.resize(n_dofs);


        // Now loop over quadrature points to handle numerical integration
        for(unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Now build the element matrix and right-hand-side using loops to
          // integrate the test functions (i) against the trial functions (j).
          for(unsigned int i=0; i<phi.size(); i++)
          {
            Fe(i) += JxW[qp]*phi[i][qp] * m_forcing;

            for(unsigned int j=0; j<phi.size(); j++)
            {
              Ke(i,j) += JxW[qp]*(this->m_paramValues(0)*dphi[i][qp]*dphi[j][qp] );
            }
          }
        }


        // At this point we have completed the matrix and RHS summation. The
        // final step is to apply boundary conditions, which in this case are
        // simple Dirichlet conditions with u(0) = u(1) = 0.

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
            Ke(s,s) += penalty;
            Fe(s)   += 0*penalty;
          }
        }

        // This is a function call that is necessary when using adaptive
        // mesh refinement. See Adaptivity Example 2 for more details.
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

        // Add Ke and Fe to the global matrix and right-hand-side.
        system.matrix->add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);
      }

      return;
    }
 


}


#endif // ASSEMBLY_VISCOUS_BURGERS_H
