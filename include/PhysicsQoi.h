#ifndef PHYSICS_QOI_H
#define PHYSICS_QOI_H


#include "agnosDefines.h"



#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/fe.h"

#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

namespace AGNOS
{


  /********************************************//**
   * \brief Derive a libmesh qoi derivative class
   ***********************************************/
  template<class T_S>
  class PhysicsQoi : public libMesh::System::QOI
  {
    public:
      PhysicsQoi(
        libMesh::EquationSystems& es, 
        const std::string&        systemName
          );
      ~PhysicsQoi( );

      void qoi( const QoISet& qoi_indices );

      void setParameterValues( const T_S& parameters );

    private:
      T_S                       m_paramValues;
      libMesh::EquationSystems& m_es; 
      const std::string&        m_systemName;

  };

  template<class T_S>
    PhysicsQoi<T_S>::PhysicsQoi( 
        libMesh::EquationSystems& es, 
        const std::string&        systemName
        )
    : m_es(es), m_systemName(systemName), m_paramValues( T_S(1) )
    { }

  template<class T_S>
    PhysicsQoi<T_S>::~PhysicsQoi( ){ }

  template<class T_S>
    void PhysicsQoi<T_S>::qoi( 
        const QoISet& qoi_indices )
    {

      // reference to mesh
      const libMesh::MeshBase& mesh = m_es.get_mesh();
      const unsigned int dim = mesh.mesh_dimension();
      libMesh::LinearImplicitSystem& system =
        m_es.get_system<LinearImplicitSystem>("1D");

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
      libMesh::Number qoiValue;

      
      libMesh::Point evalPoint(0.5);
      Q[0] = system.point_value( 0, evalPoint);


      return;
    }

  template<class T_S> 
    void PhysicsQoi<T_S>::setParameterValues( const T_S& parameters) 
    { 
      m_paramValues = parameters ;
      return; 
    }
  

}


#endif // PHYSICS_QOI_H
