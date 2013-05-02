
#ifndef PHYSICS_FUNCTION_H
#define PHYSICS_FUNCTION_H

#include "PhysicsModel.h"

namespace AGNOS
{

/********************************************//**
 * \brief Base physics function class
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunction
  {

    public:

      PhysicsFunction(
          PhysicsModel<T_S,T_P>& physics
          );
      virtual ~PhysicsFunction();

      virtual void compute( 
          const T_S& paramVector,
          T_P& imageVector
          ) = 0;

    protected:
      PhysicsModel<T_S,T_P>& m_physics;

  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunction<T_S,T_P>::PhysicsFunction( 
        PhysicsModel<T_S,T_P>& physics
        )
    : m_physics(physics)
    {
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunction<T_S,T_P>::~PhysicsFunction( )
    {
    }



/********************************************//**
 * \brief Primal solution function
 *
 *
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunctionPrimal : public PhysicsFunction<T_S,T_P>
  {
    public:
      PhysicsFunctionPrimal( PhysicsModel<T_S,T_P>& physics );
      void compute( const T_S& paramVector, T_P& imageVector) ;
  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunctionPrimal<T_S,T_P>::PhysicsFunctionPrimal( 
        PhysicsModel<T_S,T_P>& physics
        )
    : PhysicsFunction<T_S,T_P>(physics)
    { }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    void PhysicsFunctionPrimal<T_S,T_P>::compute(
        const T_S& paramVector,
        T_P& imageVector
        )
    {
      imageVector = this->m_physics.solvePrimal(paramVector);
      this->m_physics.setPrimalSolution(imageVector);
    }
  

/********************************************//**
 * \brief Adjoint solution function
 *
 *
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunctionAdjoint : public PhysicsFunction<T_S,T_P>
  {
    public:
      PhysicsFunctionAdjoint( PhysicsModel<T_S,T_P>& physics );
      void compute( const T_S& paramVector, T_P& imageVector) ;
  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunctionAdjoint<T_S,T_P>::PhysicsFunctionAdjoint( 
        PhysicsModel<T_S,T_P>& physics
        )
    : PhysicsFunction<T_S,T_P>(physics)
    { }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    void PhysicsFunctionAdjoint<T_S,T_P>::compute(
        const T_S& paramVector,
        T_P& imageVector
        )
    {
      imageVector = this->m_physics.solveAdjoint(paramVector);
    }
  


/********************************************//**
 * \brief QoI solution function
 *
 *
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunctionQoi : public PhysicsFunction<T_S,T_P>
  {
    public:
      PhysicsFunctionQoi( PhysicsModel<T_S,T_P>& physics );
      void compute( const T_S& paramVector, T_P& imageVector) ;
  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunctionQoi<T_S,T_P>::PhysicsFunctionQoi( 
        PhysicsModel<T_S,T_P>& physics
        )
    : PhysicsFunction<T_S,T_P>(physics)
    { }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    void PhysicsFunctionQoi<T_S,T_P>::compute(
        const T_S& paramVector,
        T_P& imageVector
        )
    {
      imageVector = this->m_physics.evaluateQoi(paramVector);
    }



}


#endif // PHYSICS_FUNCTION_H
