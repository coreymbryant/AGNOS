

#ifndef PHYSICS_CATENARY_H
#define PHYSICS_CATENARY_H

#include "PhysicsModel.h"

namespace AGNOS
{

  /********************************************//**
   * \brief Example PhysicsModel class - catenary chain
   *
   * A simple, single dof, system useful for testing purposes
   *
   * 
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsCatenary : public PhysicsModel<T_S,T_P>
  {

  public:
    PhysicsCatenary( );
    ~PhysicsCatenary( );

    T_P solvePrimal( 
        const T_S& parameterValue  
        );

    T_P solveAdjoint( 
        const T_S& parameterValue,  
        const T_P& primalSolution    
        );

    T_P evaluateQoi( 
        const T_S& parameterValue,
        const T_P& primalSolution    
        );

    T_P estimateError( 
        const T_S& parameterValue,  
        const T_P& primalSolution,   
        const T_P& adjointSolution  
        );
  
  protected:
    double m_T;
  };




/********************************************//**
 * \brief 
 ***********************************************/
  // NOTE: function definitions actually depend on xMin=0 and xMax=1
  //        and only coefficients for 1dof finite element solution are returned
  template<class T_S, class T_P>
  PhysicsCatenary<T_S,T_P>::PhysicsCatenary( )
    :m_T(-10.0)
  {
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsCatenary<T_S,T_P>::~PhysicsCatenary( )
  {
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsCatenary<T_S,T_P>::solvePrimal( 
        const T_S& parameterValue  
        )
    {
      T_P imageValue(1);
      

      imageValue(0) =  m_T / (8. *  parameterValue(0) ) ;
      return imageValue;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsCatenary<T_S,T_P>::solveAdjoint( 
        const T_S& parameterValue,  
        const T_P& primalSolution    
        )
    {

      T_P imageValue(1);

      imageValue(0) =  1.0 / (4. * parameterValue(0))  ;
      
      return imageValue;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsCatenary<T_S,T_P>::evaluateQoi( 
        const T_S& parameterValue,
        const T_P& primalSolution    
        )
    {
      // QoI is just primal at x=1/2 so its just coefficient value
      //return m_primalSolution ; 
      return primalSolution;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PhysicsCatenary<T_S,T_P>::estimateError( 
        const T_S& parameterValue,  
        const T_P& primalSolution,   
        const T_P& adjointSolution  
        )
    {
      // in this case the FE solution interpolates at x=1/2 so the QoI is
      // evaluated exactly
      T_P imageValue(1);
      imageValue.zero();

      return imageValue;
    }



}

#endif // PHYSICS_CATENARY_H
