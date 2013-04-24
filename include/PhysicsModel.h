
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
      ~PhysicsModel( );  /**< Default destructor */

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

      const T_P& getPrimalSolution( ) const;
      const T_P& getAdjointSolution( ) const;
      // this may need a different type
      const T_P& getErrorIndicators( ) const;

      
    protected:

      T_P* m_primalSolution;
      T_P* m_adjointSolution;
      // this may need a different type
      T_P* m_errorIndicators;
      
  }; // PhysicsModel class



  template<class T_S, class T_P>
  PhysicsModel<T_S,T_P>::PhysicsModel( )
  {
    m_primalSolution = new T_P ;
    m_adjointSolution = new T_P ;
    m_errorIndicators = new T_P ;
    return;
  }

  template<class T_S, class T_P>
  PhysicsModel<T_S,T_P>::~PhysicsModel( )
  {
    return;
  }


  template<class T_S, class T_P>
  const T_P& PhysicsModel<T_S,T_P>::getPrimalSolution( ) const
  {
    return m_primalSolution;
  }

  template<class T_S, class T_P>
  const T_P& PhysicsModel<T_S,T_P>::getAdjointSolution( ) const
  {
    return m_adjointSolution;
  }

  template<class T_S, class T_P>
  const T_P& PhysicsModel<T_S,T_P>::getErrorIndicators( ) const
  {
    return m_errorIndicators;
  }



}

#endif //PHYSICS_MODEL_H
