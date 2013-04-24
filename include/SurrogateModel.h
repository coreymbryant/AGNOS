
// should be able to build up multi-element approaches based on basic surrogate
// definition I think. Even if we use libmesh?

#ifndef SURROGATE_MODEL_H 
#define SURROGATE_MODEL_H 

namespace AGNOS
{

  /********************************************//**
   * \brief Abstract definition of the surrogate model class
   ***********************************************/
  template<class T_S, class T_P>
  class SurrogateModel
  {

    public: 

      SurrogateModel( );           /**< Default constructor */
      ~SurrogateModel( );  /**< Default destructor */

      virtual void build( ); 

      virtual T_P& evaluate( 
          T_S& parameterValues /**< parameter values to evaluate*/
          );

      virtual void refine( );

      void setParameterDimension(unsigned int parameterDimension) ;
      unsigned int getParameterDimension( ) const;

    protected: 

      PhysicsModel<T_S,T_P>* m_physicsModel ; 

      // may want to define a new class that encapsulates the paramters
      // similar to how pmpack does it
      unsigned int m_parameterDimension ;
      std::vector<double> m_parameterMins ;
      std::vector<double> m_parameterMaxs ;
      

  }; //SurrogateModel class


  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( )
    {
      return;
    }

  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::~SurrogateModel( )
    {
      return;
    }

  template<class T_S, class T_P>
    void SurrogateModel<T_S,T_P>::setParameterDimension( unsigned int parameterDimension )
    {
      m_parameterDimension = parameterDimension ;
      return;
    }

  template<class T_S, class T_P>
    unsigned int SurrogateModel<T_S,T_P>::getParameterDimension( ) const
    {
      return m_parameterDimension ;
    }

}
#endif //SURROGATE_MODEL_H
