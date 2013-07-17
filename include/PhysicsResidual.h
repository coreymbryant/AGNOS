#ifndef PHYSICS_RESIDUAL_H
#define PHYSICS_RESIDUAL_H

#include "agnosDefines.h"
#include "libmesh/nonlinear_implicit_system.h"

namespace AGNOS{

  /********************************************//**
   * \brief Base class for defining residual objects. Derived from libMesh
   * ComputeResidual class.
   ***********************************************/
  template<class T_S>
    class PhysicsResidual : public NonlinearImplicitSystem::ComputeResidual
  {

    public:
      PhysicsResidual( );
      virtual ~PhysicsResidual( );

      virtual void residual( 
          const NumericVector<Number> &X, 
          NumericVector<Number> &R,
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
    PhysicsResidual<T_S>::PhysicsResidual( )
    { }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S> 
    PhysicsResidual<T_S>::~PhysicsResidual( )
    { }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S> 
    void PhysicsResidual<T_S>::setParameterValues( const T_S& parameters) 
    { 
      m_paramValues = parameters ;
      return; 
    }





}
#endif // PHYSICS_RESIDUAL_H


