
// should be able to build up multi-element approaches based on basic surrogate
// definition I think. Even if we use libmesh?

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

      SurrogateModel( );           /**< Default constructor */
      ~SurrogateModel( );  /**< Default destructor */

      virtual void build( ); 

      virtual PhysicsDataType evaluate( 
          ParameterDataType& parameterValues /**< parameter values to evaluate*/
          );

      virtual void refine( );

      void setParameterDimension(unsigned int parameterDimension) ;
      unsigned int getParameterDimension( ) const;

    protected: 

      PhysicsModel* m_physicsModel ; 

      // may want to define a new class that encapsulates the paramters
      // similar to how pmpack does it
      unsigned int m_parameterDimension ;
      std::vector<double> m_parameterMins ;
      std::vector<double> m_parameterMaxs ;
      

  } //SurrogateModel class

}
#endif //SURROGATE_MODEL_H
