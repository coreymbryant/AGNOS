

#ifndef CATENARY_MODEL_H
#define CATENARY_MODEL_H

#define ParameterDataType double
#define PhysicsDataType double

#include "PhysicsModel.h"

namespace AGNOS
{

  template<ParameterDataType, PhysicsDataType>
  class CatenaryModel : public PhysicsModel
  {

  public:
    CatenaryModel( );
    ~CatenaryModel( );
  
  private:
    double m_T;
  };




  // NOTE: function definitions actually depend on xMin=0 and xMax=1
  //        and only coefficients for 1dof finite element solution are returned
  template<ParameterDataType, PhysicsDataType>
  CatenaryModel::CatenaryModel( )
    :m_T = -10.;
  {
    return;
  }

  template<ParameterDataType, PhysicsDataType>
    void CatenaryModel::solvePrimal( 
        ParameterDataType& parameterValue  
        )
    {
      m_primalSolution =  m_T / (8. * parameterValue) ;
      return;
    }

  template<ParameterDataType, PhysicsDataType>
    void CatenaryModel::solveAdjoint( 
        ParameterDataType& parameterValue,  
        PhysicsDataType& primalSolution    
        )
    {
      m_adjointSolution = 1.0 / (4. * parameterValue) ;
      return ;
    }

  template<ParameterDataType, PhysicsDataType>
    double CatenaryModel::evaluateQoi( 
        ParameterDataType& parameterValue,
        PhysicsDataType& primalSolution    
        )
    {
      // QoI is just primal at x=1/2 so its just coefficient value
      return m_primalSolution ; 
    }

  template<ParameterDataType, PhysicsDataType>
    double CatenaryModel::estimateError( 
        ParameterDataType& parameterValue,  
        PhysicsDataType& primalSolution,   
        PhysicsDataType& adjointSolution  
        )
    {
      // in this case the FE solution interpolates at x=1/2 so the QoI is
      // evaluated exactly
      return 0;
    }



}

#endif // CATENARY_MODEL_H
