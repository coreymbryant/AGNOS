
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
   * \brief Base surrogate model class
   *
   * Allows for derivation of Collocation and Pseudospectral surrogate models.
   * Function to construct SurrogateModel for must me defined by providing a
   * PhysicsFunction object. 
   ***********************************************/
  template<class T_S, class T_P>
  class SurrogateModel
  {


    public: 

      SurrogateModel(
          std::vector< PhysicsFunction<T_S,T_P>* >  solutionFunction,
          std::vector<Parameter*>                   parameters
          );
      SurrogateModel(
          PhysicsFunction<T_S,T_P>* solutionFunction,
          std::vector<Parameter*>   parameters
          );

      SurrogateModel( );           /**< Default constructor */
      virtual ~SurrogateModel( );  /**< Default destructor */

      // surrogate construction and evaluation
      virtual void build( ) = 0; 
      virtual std::vector<T_P> evaluate( 
          T_S& parameterValues /**< parameter values to evaluate*/
          ) = 0;
      virtual void refine( ) = 0;


      // Manipulators
      void setParameters( std::vector<Parameter*> parameters );
      std::vector<Parameter*> getParameters( ) const;

    protected: 
      
      std::vector< PhysicsFunction<T_S,T_P>* >  m_solutionFunction;
      std::vector<Parameter*>                   m_parameters;
      unsigned int                              m_dimension;

      

  }; //SurrogateModel class

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        std::vector< PhysicsFunction<T_S,T_P>* >  solutionFunction,
        std::vector<Parameter*> parameters
        )
      : m_solutionFunction(solutionFunction), m_parameters(parameters),
      m_dimension( parameters.size() )
    {
    }

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        PhysicsFunction<T_S,T_P>* solutionFunction,
        std::vector<Parameter*> parameters
        ) : m_parameters(parameters), m_dimension( parameters.size() )
    {
      m_solutionFunction = 
        std::vector< PhysicsFunction<T_S,T_P>* >(1,solutionFunction) ;
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
