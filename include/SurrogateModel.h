
// should be able to build up multi-element approaches based on basic surrogate
// definition I think. Even if we use libmesh?

#ifndef SURROGATE_MODEL_H 
#define SURROGATE_MODEL_H 

namespace AGNOS
{

  /********************************************//**
   * \brief Abstract definition of the surrogate model class
   ***********************************************/
  template<ParameterDataType, PhysicsDataType>
  class SurrogateModel
  {

    public: 

      SurrogateModel( );           /**< Default constructor */
      ~SurrogateModel( );  /**< Default destructor */

      virtual void build( ); 

      virtual PhysicsDataType& evaluate( 
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


  template<ParameterDataType, PhysicsDataType>
    SurrogateModel::SurrogateModel( )
    {
      return;
    }

  template<ParameterDataType, PhysicsDataType>
    SurrogateModel::~SurrogateModel( )
    {
      return;
    }

  template<ParameterDataType, PhysicsDataType>
    void SurrogateModel::setParameterDimension( unsigned int parameterDimension )
    {
      m_parameterDimension = parameterDimension ;
      return;
    }

  template<ParameterDataType, PhysicsDataType>
    unsigned int SurrogateModel::getParameterDimension( ) const
    {
      return m_parameterDimension ;
    }

}
#endif //SURROGATE_MODEL_H
