
#ifndef PHYSICS_MODEL_H
#define PHYSICS_MODEL_H
/********************************************//**
 * \brief Physics model class
 *
 * this is a complete description of the class
 ***********************************************/
class PhysicsModel
{
  
  public: 

    /** solve forward problem */
    virtual PhysicsDataType solvePrimal( 
        ParameterDataType parameterValue  /**< paramter value to solve at*/
        );

    /** solve adjoint problem */
    virtual PhysicsDataType solveAdjoint( 
        PhysicsDataType primalSolution,   /**< solution */
        ParameterDataType parameterValue  /**< paramter value to solve at*/
        );

    /** evaluate quantity of interest */
    virtual RealNumberDataType evaluateQoi( 
        PhysicsDataType primalSolution,   /**< solution */
        ParameterDataType parameterValue  /**< paramter value to solve at*/
        );

    /** compute physical discretization error estimate */
    virtual void estimateError( 
        PhysicsDataType primalSolution,   /**< solution */
        PhysicsDataType adjointSolution,  /**< solution */
        ParameterDataType parameterValue  /**< paramter value to solve at*/
        );


    
}

#endif //PHYSICS_MODEL_H
