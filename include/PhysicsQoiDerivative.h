#ifndef PHYSICS_QOI_DERIVATIVE_H
#define PHYSICS_QOI_DERIVATIVE_H


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
   * \brief Base QoiDerivative class for defining derivatives of Qoi for RHS of
   * adjoint problem. Derived from libmesh QOIDerivative
   ***********************************************/
  template<class T_S>
  class PhysicsQoiDerivative : public libMesh::System::QOIDerivative
  {
    public:
      PhysicsQoiDerivative(
        libMesh::EquationSystems& es, 
        const std::string&        systemName
          );
      ~PhysicsQoiDerivative( );

      virtual void qoi_derivative( const QoISet& qoi_indices ) = 0;

      void setParameterValues( const T_S& parameters );

    protected:
      T_S                       m_paramValues;
      libMesh::EquationSystems& m_es; 
      const std::string&        m_systemName;

  };

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    PhysicsQoiDerivative<T_S>::PhysicsQoiDerivative( 
        libMesh::EquationSystems& es, 
        const std::string&        systemName
        )
    : m_es(es), m_systemName(systemName), m_paramValues( T_S(1) )
    { }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S>
    PhysicsQoiDerivative<T_S>::~PhysicsQoiDerivative( ){ }


  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S> 
    void PhysicsQoiDerivative<T_S>::setParameterValues( const T_S& parameters) 
    { 
      m_paramValues = parameters ;
      return; 
    }
  

}


#endif // PHYSICS_QOI_DERIVATIVE_H
