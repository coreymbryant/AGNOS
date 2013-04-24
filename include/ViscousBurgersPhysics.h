
#ifndef VISCOUS_BURGERS_H
#define VISCOUS_BURGERS_H

#include "PhysicsModel.h"

namespace AGNOS
{
  /********************************************//**
   * \brief Basic 1D Burger's equation
   *
   * This example is given in the Le Maitre book
   ***********************************************/
  template<T_S, T_P>
  class ViscousBurgersPhysics : public PhysicsModel<T_S,T_P>
  {
    public:
      ViscousBurgersPhysics( ) ;
      ViscousBurgersPhysics(int xMin, int xMin) ;
      
      ~ViscousBurgersPhysics( ) ;

    private:
      int xMin;
      int xMax;
      double m_viscosity;
  };


  template<T_S, T_P>
    ViscousBurgersPhysics::ViscousBurgersPhysics( )
      : xMin = -10, xMax = 10
    {
      return;
    }

  template<T_S, T_P>
    ViscousBurgersPhysics::~ViscousBurgersPhysics( )
    {
      return;
    }


  template<T_S, T_P>
    void ViscousBurgersPhysics::solvePrimal(
        T_S& paramaterValue)
    {
      m_primalSolution = 1.0;
      return ;
    }


  template<T_S, T_P>
    void ViscousBurgersPhysics::solveAdjoint(
        T_S& paramaterValue,
        T_P& primalSolution
        )
    {
      m_adjointSolution = -1.0;
      return ;
    }



  template<T_S, T_P>
    double ViscousBurgersPhysics::evaluateQoi( 
        T_S& parameterValue,
        T_P& primalSolution    
        )
    {
      return 10.0;
    }

  template<T_S, T_P>
    double ViscousBurgersPhysics::estimateError( 
        T_S& parameterValue,  
        T_P& primalSolution,   
        T_P& adjointSolution  
        )
    {
      return 20.0;
    }

}

#endif // VISCOUS_BURGERS_H
