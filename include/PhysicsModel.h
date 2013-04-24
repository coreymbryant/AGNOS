
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
  template<ParameterDataType, PhysicsDataType>
  class PhysicsModel
  {
    
    public: 

      PhysicsModel( );   /**< Default constructor */
      ~PhysicsModel( );  /**< Default destructor */

      virtual void solvePrimal( 
          ParameterDataType& parameterValue  
          );

      virtual void solveAdjoint( 
          ParameterDataType& parameterValue,  
          PhysicsDataType& primalSolution    
          );

      virtual double evaluateQoi( 
          ParameterDataType& parameterValue,
          PhysicsDataType& primalSolution    
          );

      virtual double estimateError( 
          ParameterDataType& parameterValue,  
          PhysicsDataType& primalSolution,   
          PhysicsDataType& adjointSolution  
          );

      const PhysicsDataType& getPrimalSolution( ) const;
      const PhysicsDataType& getAdjointSolution( ) const;
      // this may need a different type
      const PhysicsDataType& getErrorIndicators( ) const;

      
    protected:

      PhysicsDataType& m_primalSolution;
      PhysicsDataType& m_adjointSolution;
      // this may need a different type
      PhysicsDataType& m_errorIndicators;
      
  } // PhysicsModel class

  PhysicsModel::PhysicsModel( )
  {
    return;
  }

  PhysicsModel::~PhysicsModel( )
  {
    return;
  }


  const PhysicsDataType& PhysicsModel::getPrimalSolution( ) const
  {
    return m_primalSolution;
  }

  const PhysicsDataType& PhysicsModel::getAdjointSolution( ) const
  {
    return m_adjointSolution;
  }

  const PhysicsDataType& PhysicsModel::getErrorIndicators( ) const
  {
    return m_errorIndicators;
  }



}

#endif //PHYSICS_MODEL_H
