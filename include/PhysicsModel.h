
// Do we need both versions of the defined functions, or just get by with the
// one given the solution and pass it in if needed


#ifndef PHYSICS_MODEL_H
#define PHYSICS_MODEL_H
namespace AGNOS
{

  /********************************************//**
   * \brief Abstract physics model class
   ***********************************************/
  class PhysicsModel
  {
    
    public: 

      virtual void solvePrimal( 
          ParameterDataType parameterValue  /**< paramter value to solve at*/
          );

      /** One version uses existing solution (or solves) other solves adjoint
       * for given solution */
      virtual void solveAdjoint( 
          ParameterDataType parameterValue  
          );
      virtual void solveAdjoint( 
          ParameterDataType parameterValue,  
          PhysicsDataType primalSolution    /**< solution */
          );

      /** One version uses existing solution (or solves) other evaluates QoI for
       * given solution */
      virtual RealNumberDataType evaluateQoi( 
          ParameterDataType parameterValue 
          );
      virtual RealNumberDataType evaluateQoi( 
          ParameterDataType parameterValue,
          PhysicsDataType primalSolution    /**< provided solution */
          );

      /** One version uses existing solution (or solves) other evaluates for
       * given solution */
      virtual void estimateError( 
          ParameterDataType parameterValue,  
          );
      virtual void estimateError( 
          ParameterDataType parameterValue,  
          PhysicsDataType primalSolution,   
          PhysicsDataType adjointSolution  
          );

      
    protected:
      
  } // PhysicsModel class

}

#endif //PHYSICS_MODEL_H
