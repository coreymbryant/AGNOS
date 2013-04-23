
#ifndef SURROGATE_MODEL_H 
#define SURROGATE_MODEL_H 

namespace AGNOS
{

  /********************************************//**
   * \brief Abstract definition of the surrogate model class
   ***********************************************/
  class SurrogateModel
  {

    public: 

      SurrogateModel();           /**< Default constructor */
      virtual ~SurrogateModel();  /**< Default destructor */

      virtual void build( ); 

      virtual PhysicsDataType evaluate( 
          ParameterDataType& parameterValues /**< parameter values to evaluate*/
          );

      virtual void refine( );

      void setPhysicsModel(PhysicsModel& physicsModel) const;
      const PhysicsModel& getPhysicsModel( ) const;

    protected: 

      const PhysicsModel& m_physicsModel ; 

  } //SurrogateModel class

}
#endif //SURROGATE_MODEL_H
