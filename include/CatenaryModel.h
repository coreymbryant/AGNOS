

#ifndef CATENARY_MODEL_H
#define CATENARY_MODEL_H

#include "PhysicsModel.h"

namespace AGNOS
{

  template<class T_S, class T_P>
  class CatenaryModel : public PhysicsModel<T_S,T_P>
  {

  public:
    CatenaryModel( );
    ~CatenaryModel( );

    virtual void solvePrimal( 
        T_S& parameterValue  
        );

    virtual void solveAdjoint( 
        T_S& parameterValue,  
        T_P& primalSolution    
        );

    virtual double evaluateQoi( 
        T_S& parameterValue,
        T_P& primalSolution    
        );

    virtual double estimateError( 
        T_S& parameterValue,  
        T_P& primalSolution,   
        T_P& adjointSolution  
        );
  
  private:
    double m_T;
  };




  // NOTE: function definitions actually depend on xMin=0 and xMax=1
  //        and only coefficients for 1dof finite element solution are returned
  template<class T_S, class T_P>
  CatenaryModel<T_S,T_P>::CatenaryModel( )
    :m_T(-10.0)
  {
    return;
  }

  template<class T_S, class T_P>
    void CatenaryModel<T_S,T_P>::solvePrimal( 
        T_S& parameterValue  
        )
    {
      //m_primalSolution =  m_T / (8. * parameterValue) ;
      return;
    }

  template<class T_S, class T_P>
    void CatenaryModel<T_S,T_P>::solveAdjoint( 
        T_S& parameterValue,  
        T_P& primalSolution    
        )
    {
      //m_adjointSolution = 1.0 / (4. * parameterValue) ;
      return ;
    }

  template<class T_S, class T_P>
    double CatenaryModel<T_S,T_P>::evaluateQoi( 
        T_S& parameterValue,
        T_P& primalSolution    
        )
    {
      // QoI is just primal at x=1/2 so its just coefficient value
      //return m_primalSolution ; 
      return 0.;
    }

  template<class T_S, class T_P>
    double CatenaryModel<T_S,T_P>::estimateError( 
        T_S& parameterValue,  
        T_P& primalSolution,   
        T_P& adjointSolution  
        )
    {
      // in this case the FE solution interpolates at x=1/2 so the QoI is
      // evaluated exactly
      return 0;
    }



}

#endif // CATENARY_MODEL_H
