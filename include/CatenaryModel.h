

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
  
  private:
    double m_T;
  };




/********************************************//**
 * \brief 
 ***********************************************/
  // NOTE: function definitions actually depend on xMin=0 and xMax=1
  //        and only coefficients for 1dof finite element solution are returned
  template<class T_S, class T_P>
  CatenaryModel<T_S,T_P>::CatenaryModel( )
    :m_T(-10.0)
  {
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  CatenaryModel<T_S,T_P>::~CatenaryModel( )
  {
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P CatenaryModel<T_S,T_P>::solvePrimal( 
        const T_S& parameterValue  
        )
    {
      T_P imageValue(1,0.0);

      imageValue[0] =  m_T / (8. * parameterValue[0] ) ;
      *(this->m_primalSolution) = imageValue ;
      return imageValue;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P CatenaryModel<T_S,T_P>::solveAdjoint( 
        const T_S& parameterValue,  
        const T_P& primalSolution    
        )
    {

      T_P imageValue(1,0.0);

      imageValue[0] =  1.0 / (4. * parameterValue[0])  ;
      *(this->m_adjointSolution) = imageValue ;
      
      return imageValue;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P CatenaryModel<T_S,T_P>::evaluateQoi( 
        const T_S& parameterValue,
        const T_P& primalSolution    
        )
    {
      // QoI is just primal at x=1/2 so its just coefficient value
      //return m_primalSolution ; 
      T_P dummy;
      return dummy;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    T_P CatenaryModel<T_S,T_P>::estimateError( 
        const T_S& parameterValue,  
        const T_P& primalSolution,   
        const T_P& adjointSolution  
        )
    {
      // in this case the FE solution interpolates at x=1/2 so the QoI is
      // evaluated exactly
      T_P dummy;
      return dummy;
    }



}

#endif // CATENARY_MODEL_H
