
// Do we need both versions of the defined functions, or just get by with the
// one given the solution and pass it in if needed
//
// Solve functions intentially do not return the solution. There may be cases
// where we want to solve the problem internally but not create a new instance
// of the solution. Create get..Solution functions to check if solution exists
// and if not solve problem first. These should be what is used by surrogate
// model.



#ifndef PHYSICS_MODEL_H
#define PHYSICS_MODEL_H

#include "libmesh/error_vector.h"

namespace AGNOS
{

  /********************************************//**
   * \brief Base physics model class
   *
   * Abstract framework for PhysicsModel classes. 
   ***********************************************/
  template<class T_S, class T_P>
  class PhysicsModel
  {
    
    public: 

      PhysicsModel( );   /**< Default constructor */
      virtual ~PhysicsModel( );  /**< Default destructor */

      virtual T_P solvePrimal( 
          const T_S& parameterValue  
          ) = 0;
      void setPrimalSolution( T_P primalSolution );

      virtual T_P solveAdjoint( 
          const T_S& parameterValue,  
          const T_P& primalSolution    
          ) = 0;
      void setAdjointSolution( T_P adjointSolution );

      virtual T_P evaluateQoi( 
          const T_S& parameterValue,
          const T_P& primalSolution    
          ) = 0;
      void setQoiValue( T_P qoiValue );

      virtual T_P estimateError( 
          const T_S& parameterValue,  
          const T_P& primalSolution,   
          const T_P& adjointSolution  
          ) = 0;
      void setErrorEstimate( T_P errorEstimate );
      void setErrorIndicators( T_P errorIndicators );

      virtual void refine (  )  = 0;
      virtual void refine ( T_P& errorIndicators  ) ;

      const T_P* getPrimalSolution( ) const ;
      const T_P* getAdjointSolution( ) const ;
      const T_P* getQoiValue( ) const;
      const T_P* getErrorEstimate( ) const;
      const T_P* getErrorIndicators( ) const;

      void resetSolution( );

      bool useUniformRefinement();
      bool resolveAdjoint();
      
    protected:

      T_P* m_primalSolution;
      T_P* m_adjointSolution;
      T_P* m_qoiValue;
      T_P* m_errorEstimate;
      T_P* m_errorIndicators;
      T_P* m_meanErrorIndicator;

      bool m_useUniformRefinement;
      
      // for nonlinear problems we need to resolve the adjoint at each
      // parameterValue in TotalErrorFunction 
      bool m_resolveAdjoint;

      
  }; // PhysicsModel class



/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsModel<T_S,T_P>::PhysicsModel( )
  : m_resolveAdjoint(false),m_useUniformRefinement(true)
  {
    m_primalSolution  = NULL;
    m_adjointSolution = NULL;
    m_qoiValue        = NULL;
    m_errorEstimate = NULL;
    m_errorIndicators = NULL;

  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsModel<T_S,T_P>::~PhysicsModel( )
  {
    if ( m_primalSolution )
      delete m_primalSolution;
    if ( m_adjointSolution )
      delete m_adjointSolution;
    if ( m_qoiValue )
      delete m_qoiValue;
    if ( m_errorEstimate )
      delete m_errorEstimate;
    if ( m_errorIndicators )
      delete m_errorIndicators;
  }


/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  const T_P* PhysicsModel<T_S,T_P>::getPrimalSolution() const
  {
    return m_primalSolution;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  const T_P* PhysicsModel<T_S,T_P>::getAdjointSolution( ) const
  {
    return m_adjointSolution;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  const T_P* PhysicsModel<T_S,T_P>::getQoiValue( ) const
  {
    return m_qoiValue;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  const T_P* PhysicsModel<T_S,T_P>::getErrorEstimate( ) const
  {
    return m_errorEstimate;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  const T_P* PhysicsModel<T_S,T_P>::getErrorIndicators( ) const
  {
    return m_errorIndicators;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsModel<T_S,T_P>::setPrimalSolution( 
      T_P primalSolution )
  {
    delete m_primalSolution;
    m_primalSolution = new T_P(primalSolution);
    return ;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsModel<T_S,T_P>::setAdjointSolution( 
      T_P adjointSolution )
  {
    delete m_adjointSolution;
    m_adjointSolution = new T_P(adjointSolution);
    return ;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsModel<T_S,T_P>::setQoiValue( 
      T_P qoiValue )
  {
    delete m_qoiValue;
    m_qoiValue = new T_P(qoiValue);
    return ;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsModel<T_S,T_P>::setErrorEstimate( 
      T_P errorEstimate )
  {
    delete m_errorEstimate;
    m_errorEstimate = new T_P(errorEstimate);
    return ;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsModel<T_S,T_P>::setErrorIndicators( 
      T_P errorIndicators )
  {
    delete m_errorIndicators;
    m_errorIndicators = new T_P(errorIndicators);
    return ;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsModel<T_S,T_P>::resetSolution( )
  {
    if ( m_primalSolution )
      delete m_primalSolution;
    if ( m_adjointSolution )
      delete m_adjointSolution;
    if ( m_qoiValue )
      delete m_qoiValue;
    if ( m_errorEstimate )
      delete m_errorEstimate;
    if ( m_errorIndicators )
      delete m_errorIndicators;

    m_primalSolution  = NULL;
    m_adjointSolution = NULL;
    m_qoiValue        = NULL;
    m_errorEstimate = NULL;
    m_errorIndicators = NULL;
  }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
    bool PhysicsModel<T_S,T_P>::useUniformRefinement()
    {
      return m_useUniformRefinement;
    }


  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
    void PhysicsModel<T_S,T_P>::refine(
        T_P& errorIndicators)
    {
      std::cerr 
        << "\n\t ERROR: Adaptive refinement is not available in derived "
        << "PhysicsModel class\n" << std::endl;
      exit(1);
      return;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  bool PhysicsModel<T_S,T_P>::resolveAdjoint()
  {
    return m_resolveAdjoint;
  }

}

#endif //PHYSICS_MODEL_H
