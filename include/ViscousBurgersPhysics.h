
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
  template<ParameterDataType, PhysicsDataType>
  class ViscousBurgersPhysics : public PhysicsModel
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


  template<ParameterDataType, PhysicsDataType>
    ViscousBurgersPhysics::ViscousBurgersPhysics( )
      : xMin = -10, xMax = 10
    {
      return;
    }

  template<ParameterDataType, PhysicsDataType>
    ViscousBurgersPhysics::~ViscousBurgersPhysics( )
    {
      return;
    }


  template<ParameterDataType, PhysicsDataType>
    void ViscousBurgersPhysics::solvePrimal(
        ParameterDataType& paramaterValue)
    {
      m_primalSolution = 1.0;
      return ;
    }


  template<ParameterDataType, PhysicsDataType>
    void ViscousBurgersPhysics::solveAdjoint(
        ParameterDataType& paramaterValue,
        PhysicsDataType& primalSolution
        )
    {
      m_adjointSolution = -1.0;
      return ;
    }



  template<ParameterDataType, PhysicsDataType>
    double ViscousBurgersPhysics::evaluateQoi( 
        ParameterDataType& parameterValue,
        PhysicsDataType& primalSolution    
        )
    {
      return 10.0;
    }

  template<ParameterDataType, PhysicsDataType>
    double ViscousBurgersPhysics::estimateError( 
        ParameterDataType& parameterValue,  
        PhysicsDataType& primalSolution,   
        PhysicsDataType& adjointSolution  
        )
    {
      return 20.0;
    }

}

#endif // VISCOUS_BURGERS_H
