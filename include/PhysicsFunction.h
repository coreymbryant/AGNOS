
#ifndef PHYSICS_FUNCTION_H
#define PHYSICS_FUNCTION_H

#include <map>

#include "PhysicsModel.h"

namespace AGNOS
{

/********************************************//**
 * \brief Base physics function class
 *
 * This class provides a means for accessing primal and adjoint solutions from
 * the PhysicsModel class. It also allows the definition of simple
 * (user-defined) functions that are independent of a PhysicsModel object. 
 *
 * It provides a convenient, uniform, interface for a SurrogateModel object to
 * call, namely the compute() function. 
 * 
 * 
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunction
  {

    public:

      PhysicsFunction();
      virtual ~PhysicsFunction();

      virtual void compute( 
          const T_S& paramVector,
          T_P& imageVector
          ) = 0;

      std::string name( ) const;

    protected:
      std::string m_name;

  };



/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunction<T_S,T_P>::PhysicsFunction( ) 
    { }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunction<T_S,T_P>::~PhysicsFunction( )
    { }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P> 
    std::string PhysicsFunction<T_S,T_P>::name( ) const
    {
      return m_name;
    }
  



/********************************************//**
 * \brief Basic function independent of a PhysicsModel object
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunctionSimple : public PhysicsFunction<T_S,T_P>
  {
    public:
      PhysicsFunctionSimple( 
          std::string functionName,
          T_P (*myFunction)(const T_S&)  
          ) : m_myFunction(myFunction) 
      { this->m_name = functionName; }

      void compute( const T_S& paramVector, T_P& imageVector)
      {
        imageVector = m_myFunction(paramVector) ;
        return ;
      };

    protected:
      T_P (*m_myFunction)(const T_S&);
  };





/********************************************//**
 * \brief Primal solution function based on a PhysicsModel object
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunctionPrimal : public PhysicsFunction<T_S,T_P>
  {
    public:
      PhysicsFunctionPrimal( PhysicsModel<T_S,T_P>& physics );
      void compute( const T_S& paramVector, T_P& imageVector) ;
    protected:
      PhysicsModel<T_S,T_P>& m_physics;
  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunctionPrimal<T_S,T_P>::PhysicsFunctionPrimal( 
        PhysicsModel<T_S,T_P>& physics
        )
    : m_physics(physics)
    { this->m_name = "primal"; }

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
      return ;
    }
  

/********************************************//**
 * \brief Adjoint solution function based on a PhysicsModel object
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunctionAdjoint : public PhysicsFunction<T_S,T_P>
  {
    public:
      PhysicsFunctionAdjoint( PhysicsModel<T_S,T_P>& physics );
      void compute( const T_S& paramVector, T_P& imageVector) ;
    protected:
      PhysicsModel<T_S,T_P>& m_physics;
  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunctionAdjoint<T_S,T_P>::PhysicsFunctionAdjoint( 
        PhysicsModel<T_S,T_P>& physics
        )
    : m_physics(physics)
    { this->m_name = "adjoint"; }

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
 * \brief QoI solution function based on a PhysicsModel object
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunctionQoi : public PhysicsFunction<T_S,T_P>
  {
    public:
      PhysicsFunctionQoi( PhysicsModel<T_S,T_P>& physics );
      void compute( const T_S& paramVector, T_P& imageVector) ;
    protected:
      PhysicsModel<T_S,T_P>& m_physics;
  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunctionQoi<T_S,T_P>::PhysicsFunctionQoi( 
        PhysicsModel<T_S,T_P>& physics
        )
    : m_physics(physics)
    { this->m_name = "qoi"; }

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
