
// should be able to build up multi-element approaches based on basic surrogate
// definition I think. Even if we use libmesh?

#ifndef SURROGATE_MODEL_H 
#define SURROGATE_MODEL_H 
#include "Parameter.h"
#include "PhysicsFunction.h"
#include <iostream>

namespace AGNOS
{

  /********************************************//**
   * \brief Abstract definition of the surrogate model class
   ***********************************************/
  template<class T_S, class T_P>
  class SurrogateModel
  {


    public: 

      // provided physics must me defined
      // easiest way to accomplish this is to pass a function pointer
      // alternatively we could have made a virtual function, but we would have
      // to derive many additional classes to group forward and adjoint solution
      // evaluations etc. 
      SurrogateModel(
          PhysicsFunction<T_S,T_P>& solutionFunction,
          std::vector<Parameter*>   parameters
          );

      SurrogateModel( );           /**< Default constructor */
      virtual ~SurrogateModel( );  /**< Default destructor */

      // surrogate construction and evaluation
      virtual void build( ) = 0; 
      virtual T_P evaluate( 
          T_S& parameterValues /**< parameter values to evaluate*/
          ) = 0;
      virtual void refine( ) = 0;


      // Manipulators
      void setParameters( std::vector<Parameter*> parameters );
      std::vector<Parameter*> getParameters( ) const;

    protected: 
      
      PhysicsFunction<T_S,T_P>& m_solutionFunction;
      std::vector<Parameter*>   m_parameters;
      unsigned int              m_dimension;

      

  }; //SurrogateModel class

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        PhysicsFunction<T_S,T_P>& solutionFunction,
        std::vector<Parameter*> parameters
        )
      : m_solutionFunction(solutionFunction), m_parameters(parameters),
      m_dimension( parameters.size() )
    {
    }

/********************************************//*
 * \brief Default constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( )
    {
    }

/********************************************//**
 * \brief Default destructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::~SurrogateModel( )
    {
    }

/********************************************//**
 * \brief set parameters
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModel<T_S,T_P>::setParameters( 
        std::vector<Parameter*> parameters)
    {
      m_parameters = parameters ;
      m_dimension = parameters.size();
      return;
    }

/********************************************//**
 * \brief get number of parameters
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<Parameter*> SurrogateModel<T_S,T_P>::getParameters( ) const
    {
      return m_parameters ;
    }

}
#endif //SURROGATE_MODEL_H
