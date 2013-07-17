#ifndef PHYSICS_JACOBIAN_H
#define PHYSICS_JACOBIAN_H

#include "agnosDefines.h"
#include "libmesh/nonlinear_implicit_system.h"

namespace AGNOS{

  /********************************************//**
   * \brief Base class for defining jacobian objects. Derived from libMesh
   * ComputeJacobian class.
   ***********************************************/
  template<class T_S>
    class PhysicsJacobian : public NonlinearImplicitSystem::ComputeJacobian
  {

    public:
      PhysicsJacobian( );
      virtual ~PhysicsJacobian( );

      virtual void jacobian( 
          const NumericVector<Number> &X, 
          SparseMatrix<Number> &J,
          NonlinearImplicitSystem &S) = 0;

      void setParameterValues( const T_S& parameters );
      virtual void setSystemData( const GetPot& input ) = 0;

    protected:
      T_S m_paramValues;
  
  };

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S> 
    PhysicsJacobian<T_S>::PhysicsJacobian()
    {}

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S> 
    PhysicsJacobian<T_S>::~PhysicsJacobian()
    {}

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S> 
    void PhysicsJacobian<T_S>::setParameterValues( const T_S& parameters) 
    { 
      m_paramValues = parameters ;
      return; 
    }





}
#endif // PHYSICS_JACOBIAN_H


