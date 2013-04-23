
#include "ViscousBurgersPhysics.h"

namespace AGNOS
{

  ViscousBurgersPhysics::ViscousBurgersPhysics( )
    : xMin = -10, xMax = 10
  {
    return;
  }

  ViscousBurgersPhysics::~ViscousBurgersPhysics( )
  {
    return;
  }


  void ViscousBurgersPhysics::solvePrimal(
      ParameterDataType& paramaterValue)
  {
    m_primalSolution = 1.0;
    return ;
  }


  void ViscousBurgersPhysics::solveAdjoint(
      ParameterDataType& paramaterValue,
      PhysicsDataType& primalSolution
      )
  {
    m_adjointSolution = -1.0;
    return ;
  }



    double ViscousBurgersPhysics::evaluateQoi( 
        ParameterDataType& parameterValue,
        PhysicsDataType& primalSolution    
        )
    {
      return 10.0;
    }

    double ViscousBurgersPhysics::estimateError( 
        ParameterDataType& parameterValue,  
        PhysicsDataType& primalSolution,   
        PhysicsDataType& adjointSolution  
        )
    {
      return 20.0;
    }

}
