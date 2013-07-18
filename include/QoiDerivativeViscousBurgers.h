#ifndef QOI_DERIVATIVE_VISCOUS_BURGERS_H
#define QOI_DERIVATIVE_VISCOUS_BURGERS_H


#include "agnosDefines.h"

#include "PhysicsQoiDerivative.h"

namespace AGNOS
{


  /********************************************//**
   * \brief Derive a libmesh qoi derivative class
   ***********************************************/
  template<class T_S>
  class QoiDerivativeViscousBurgers : public PhysicsQoiDerivative<T_S>
  {
    public:
      QoiDerivativeViscousBurgers(
        libMesh::EquationSystems& es, 
        const std::string&        systemName
          );
      ~QoiDerivativeViscousBurgers( );

      void qoi_derivative( const QoISet& qoi_indices );


    private:

  };

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    QoiDerivativeViscousBurgers<T_S>::QoiDerivativeViscousBurgers( 
        libMesh::EquationSystems& es, 
        const std::string&        systemName
        )
    : PhysicsQoiDerivative<T_S>(es,systemName)
    { }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    QoiDerivativeViscousBurgers<T_S>::~QoiDerivativeViscousBurgers( ){ }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    void QoiDerivativeViscousBurgers<T_S>::qoi_derivative( 
        const QoISet& qoi_indices )
    {

      // reference to mesh
      const libMesh::MeshBase& mesh = this->m_es.get_mesh();
      const unsigned int dim = mesh.mesh_dimension();
      libMesh::System& system = this->m_es.get_system("Burgers");

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
      const std::vector<Point > &q_point = fe->get_xyz();

      // declare matrices
      /* libMesh::DenseMatrix<libMesh::Number> Ke; */
      /* libMesh::DenseVector<libMesh::Number> Fe; */
      libMesh::DenseVector<libMesh::Number> Qe;
      libMesh::NumericVector<libMesh::Number> &Q = system.get_adjoint_rhs();
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
        /* Ke.resize(n_dofs, n_dofs); */
        /* Fe.resize(n_dofs); */
        Qe.resize(n_dofs);


        // Now loop over quadrature points to handle numerical integration
        for(unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
           const Real x = q_point[qp](0);

          // Now build the element matrix and right-hand-side using loops to
          // integrate the test functions (i) against the trial functions (j).
          for(unsigned int i=0; i<phi.size(); i++)
          {
            Qe(i) += JxW[qp]*phi[i][qp] ;
          }
        }

        Q.add_vector(Qe, dof_indices);

      }

      return;
    }

  

}


#endif // QOI_DERIVATIVE_VISCOUS_BURGERS_H
