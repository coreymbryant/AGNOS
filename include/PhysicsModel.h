
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
namespace AGNOS
{

  /********************************************//**
   * \brief Abstract physics model class
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

      virtual T_P solveAdjoint( 
          const T_S& parameterValue,  
          const T_P& primalSolution    
          ) = 0;
      virtual T_P solveAdjoint( const T_S& parameterValue ) ;

      virtual T_P evaluateQoi( 
          const T_S& parameterValue,
          const T_P& primalSolution    
          ) = 0;
      virtual T_P evaluateQoi( const T_S& parameterValue ) ;

      virtual T_P estimateError( 
          const T_S& parameterValue,  
          const T_P& primalSolution,   
          const T_P& adjointSolution  
          ) = 0;
      virtual T_P estimateError( const T_S& parameterValue ) ;

      const T_P* getPrimalSolution( ) ;
      const T_P* getAdjointSolution( ) ;

      // this may need a different type
      const T_P* getErrorIndicators( ) ;

      
    protected:

      T_P* m_primalSolution;
      T_P* m_adjointSolution;
      T_P* m_qoiValue;
      // this may need a different type
      T_P* m_errorIndicators;

      
  }; // PhysicsModel class



/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsModel<T_S,T_P>::PhysicsModel( )
  {
    m_primalSolution  (NULL);
    m_adjointSolution (NULL);
    m_qoiValue        (NULL);
    m_errorIndicators (NULL);
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  PhysicsModel<T_S,T_P>::~PhysicsModel( )
  {
    delete [] m_primalSolution;
    delete [] m_adjointSolution;
    delete [] m_errorIndicators;
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
  const T_P PhysicsModel<T_S,T_P>::getErrorIndicators( ) const
  {
    return m_errorIndicators;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  T_P PhysicsModel<T_S,T_P>::solveAdjoint( const T_S& parameterValue ) ;

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  virtual T_P PhysicsModel<T_S,T_P>::evaluateQoi( const T_S& parameterValue ) ;

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
  virtual T_P PhysicsModel<T_S,T_P>::estimateError( const T_S& parameterValue ) ;

}

#endif //PHYSICS_MODEL_H
