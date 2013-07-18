#ifndef QOI_CATENARY_H
#define QOI_CATENARY_H


#include "agnosDefines.h"

#include "PhysicsQoi.h"

namespace AGNOS
{


  /********************************************//**
   * \brief Derived PhysicsQoi class for catenary problem
   ***********************************************/
  template<class T_S>
  class QoiCatenary : public PhysicsQoi<T_S>
  {
    public:
      QoiCatenary(
        libMesh::EquationSystems& es, 
        const std::string&        systemName
          );
      ~QoiCatenary( );

      void qoi( const QoISet& qoi_indices );

    private:

  };

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    QoiCatenary<T_S>::QoiCatenary( 
        libMesh::EquationSystems& es, 
        const std::string&        systemName
        )
    : PhysicsQoi<T_S>( es, systemName )
    { }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    QoiCatenary<T_S>::~QoiCatenary( ){ }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    void QoiCatenary<T_S>::qoi( const QoISet& qoi_indices )
    {

      // reference to mesh
      const libMesh::MeshBase& mesh = this->m_es.get_mesh();
      const unsigned int dim = mesh.mesh_dimension();
      libMesh::System& system = this->m_es.get_system("1D");

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
      std::vector<libMesh::Number> &Q = system.qoi;
      Q[0]=0;
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



        // Now loop over quadrature points to handle numerical integration
        for(unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
           const Real x = q_point[qp](0);

           // get solution value 
           Number U = system.point_value(0,q_point[qp]);

          // Now build the element matrix and right-hand-side using loops to
          // integrate the test functions (i) against the trial functions (j).
            Q[0] += JxW[qp] * U
              * sqrt(100./(3.14159)) * std::exp( -100. * pow(x-0.5,2) );
        }


      }

      // overwrite to exact point calculation
      libMesh::Point point(0.5);
      Q[0] = system.point_value(0,point);

      return;
    }


}


#endif // QOI_CATENARY_H