
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
      void setPrimalSolution( T_P& primalSolution );

      virtual T_P solveAdjoint( 
          const T_S& parameterValue,  
          const T_P& primalSolution    
          ) = 0;
      const T_P solveAdjoint( const T_S& parameterValue) ;

      virtual T_P evaluateQoi( 
          const T_S& parameterValue,
          const T_P& primalSolution    
          ) = 0;
      const T_P evaluateQoi( const T_S& parameterValue);

      virtual T_P estimateError( 
          const T_S& parameterValue,  
          const T_P& primalSolution,   
          const T_P& adjointSolution  
          ) = 0;
      const T_P estimateError(  const T_S& parameterValue );

      virtual void refine (  )  = 0;

      const T_P* getPrimalSolution( ) const;
      const T_P* getAdjointSolution( ) const;
      const T_P* getErrorIndicators( ) const;


      
    protected:

      T_P* m_primalSolution;
      T_P* m_adjointSolution;
      T_P* m_qoiValue;
      // this may need a different type
      T_P*                            m_errorEstimate;
      T_P*                            m_errorIndicators;
      libMesh::ErrorVector            m_meanErrorIndicator;

      
  }; // PhysicsModel class



/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsModel<T_S,T_P>::PhysicsModel( )
  {
    /* m_primalSolution  = new T_P; */
    /* m_adjointSolution = new T_P; */
    /* m_qoiValue        = new T_P; */
    /* m_errorIndicators = new T_P; */

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
  const T_P* PhysicsModel<T_S,T_P>::getErrorIndicators( ) const
  {
    return m_errorIndicators;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  void PhysicsModel<T_S,T_P>::setPrimalSolution( 
      T_P& primalSolution )
  {
    delete m_primalSolution;
    m_primalSolution = new T_P(primalSolution);
    /* *m_primalSolution = primalSolution ; */

    return ;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  const T_P PhysicsModel<T_S,T_P>::solveAdjoint( 
      const T_S& parameterValue )
  {
    //TODO we could have problems if sol is for different parameter value
    if (m_primalSolution == NULL)
    {
      delete m_primalSolution;
      m_primalSolution = new T_P( solvePrimal( parameterValue ) ) ;
    }

    delete m_adjointSolution;
    m_adjointSolution = new T_P( solveAdjoint( parameterValue, *m_primalSolution) );

    return *m_adjointSolution;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  const T_P PhysicsModel<T_S,T_P>::evaluateQoi( 
      const T_S& parameterValue )
  {

    if (m_primalSolution == NULL)
    {
      delete m_primalSolution;
      m_primalSolution = new T_P( solvePrimal( parameterValue ) ) ;
    }

    delete m_qoiValue;
    m_qoiValue = new T_P( evaluateQoi( parameterValue, *m_primalSolution ) );


    return *m_qoiValue;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  const T_P PhysicsModel<T_S,T_P>::estimateError( 
      const T_S& parameterValue )
  {

    if (m_primalSolution == NULL)
    {
      delete m_primalSolution;
      m_primalSolution = new T_P( solvePrimal( parameterValue ) ) ;
      m_adjointSolution = new T_P( solveAdjoint( parameterValue ) ) ;
    }

    delete m_qoiValue;
    m_qoiValue = new T_P( evaluateQoi( parameterValue, *m_primalSolution ) );
    m_errorIndicators = new T_P( estimateError( parameterValue, *m_primalSolution,
          *m_adjointSolution ) );

    return *m_errorIndicators;
  }

}

#endif //PHYSICS_MODEL_H
