

/********************************************//**
 * Physics model class
 ***********************************************/
class PhysicsModel
{
  
  public: 

    /** solve forward problem */
    virtual PhysicsDataType solvePrimal( 
        ParameterDataType parameterValue /**< paramter value to solve at*/
        );

    
}
