
#ifndef PHYSICS_VISCOUS_BURGERS_H
#define PHYSICS_VISCOUS_BURGERS_H

#include "PhysicsModel.h"

namespace AGNOS
{
  /********************************************//**
   * \brief Basic 1D Burger's equation
   *
   * This example is given in the Le Maitre book
   ***********************************************/
  template<T_S, T_P>
  class PhysicsViscousBurgers : public PhysicsModel<T_S,T_P>
  {
    public:
      PhysicsViscousBurgers( ) ;
      PhysicsViscousBurgers(int xMin, int xMin) ;
      
      ~PhysicsViscousBurgers( ) ;

    private:
      int xMin;
      int xMax;
      double m_viscosity;
  };


  template<T_S, T_P>
    PhysicsViscousBurgers::PhysicsViscousBurgers( )
      : xMin = -10, xMax = 10
    {
      return;
    }

  template<T_S, T_P>
    PhysicsViscousBurgers::~PhysicsViscousBurgers( )
    {
      return;
    }


  template<T_S, T_P>
    void PhysicsViscousBurgers::solvePrimal(
        T_S& paramaterValue)
    {
      m_primalSolution = 1.0;
      return ;
    }


  template<T_S, T_P>
    void PhysicsViscousBurgers::solveAdjoint(
        T_S& paramaterValue,
        T_P& primalSolution
        )
    {
      m_adjointSolution = -1.0;
      return ;
    }



  template<T_S, T_P>
    double PhysicsViscousBurgers::evaluateQoi( 
        T_S& parameterValue,
        T_P& primalSolution    
        )
    {
      return 10.0;
    }

  template<T_S, T_P>
    double PhysicsViscousBurgers::estimateError( 
        T_S& parameterValue,  
        T_P& primalSolution,   
        T_P& adjointSolution  
        )
    {
      return 20.0;
    }

}

#endif // PHYSICS_VISCOUS_BURGERS_H
