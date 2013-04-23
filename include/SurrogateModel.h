
#ifndef SURROGAT_MODEL_H 
#define SURROGAT_MODEL_H 

/********************************************//**
 * Abstract definition of the surrogate model class
 ***********************************************/
class SurrogateModel
{

  public: 

    SurrogateModel();           /**< Default constructor */
    virtual ~SurrogateModel();  /**< Default destructor */

    /** construct surrogate model */
    virtual void construct( 
        PhysicsModel& physicsModel /**< underlying physics model class */
        ); 

    /** evaluate surrogate model */
    virtual PhysicsDataType evaluate( 
        ParameterDataType& parameterValues /**< parameter values to evaluate at*/
        );

    /**< refine surrogate model*/
    virtual void refine( );

};

#endif //SURROGAT_MODEL_H
