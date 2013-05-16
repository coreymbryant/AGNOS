
// should be able to build up multi-element approaches based on basic surrogate
// definition I think. Even if we use libmesh?

#ifndef SURROGATE_MODEL_H 
#define SURROGATE_MODEL_H 
#include "Parameter.h"
#include "PhysicsFunction.h"
#include <iostream>

namespace AGNOS
{

  enum SurrogateModelType{
    PSEUDO_SPECTRAL_TENSOR_PRODUCT=0,
    PSEUDO_SPECTRAL_SPARSE_GRID,
    PSEUDO_SPECTRAL_MONTE_CARLO,
    COLLOCATION };

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

      // single physics function constructors
      SurrogateModel(
          PhysicsFunction<T_S,T_P>*         solutionFunction,
          const std::vector<Parameter*>     parameters,
          const unsigned int                order 
          );
      SurrogateModel(
          PhysicsFunction<T_S,T_P>*         solutionFunction,
          const std::vector<Parameter*>     parameters,
          const std::vector<unsigned int>&  order
          );

      // multiple physics function constructors
      SurrogateModel(
          std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
          const std::vector<Parameter*>                       parameters,
          const unsigned int                                  order 
          );
      SurrogateModel(
          std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
          const std::vector<Parameter*>                       parameters,
          const std::vector<unsigned int>&                    order
          );

      SurrogateModel( );           /**< Default constructor */
      virtual ~SurrogateModel( );  /**< Default destructor */

      // surrogate construction and evaluation
      virtual void build( ) = 0; 
      virtual std::map<std::string, T_P> evaluate( 
          T_S& parameterValues /**< parameter values to evaluate*/
          ) = 0;
      virtual void refine( ) = 0;


      // Manipulators
      void setParameters( std::vector<Parameter*> parameters );
      std::vector<Parameter*> getParameters( ) const;
      const std::map< std::string, std::vector<T_P> >   
                                getCoefficients( ) const;
      std::vector<unsigned int> getExpansionOrder( ) const;

    protected: 
      
      std::vector<unsigned int>                           m_order;  
      std::map< std::string, std::vector<T_P> >           m_coefficients;

      std::map< std::string, PhysicsFunction<T_S,T_P>* >  m_solutionFunction;
      std::vector<Parameter*>                             m_parameters;
      unsigned int                                        m_dimension;

      

  }; //SurrogateModel class

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
        std::vector<Parameter*>                             parameters,
        unsigned int                                        order
        )
      : 
        m_solutionFunction(solutionFunction), 
        m_parameters(parameters),
        m_dimension( parameters.size() )
    {
      m_order = std::vector<unsigned int>(m_dimension,order);
      // TODO how to do this
      /* m_coefficients.resize( solutionFunction.size() ); */
    }

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
        std::vector<Parameter*>                   parameters,
        const std::vector<unsigned int>&          order
        )
      : 
        m_solutionFunction(solutionFunction), 
        m_parameters(parameters),
        m_dimension( parameters.size() ), 
        m_order(order)
    {
      /* m_coefficients.resize( solutionFunction.size() ); */

      if (m_order.size() != parameters.size() )
      {
        std::cout 
          << std::endl
          << "\tERROR:"
          << " order vector dimension does not match number of parameters"
          << std::endl
          << std::endl;
        assert(0);
      }
    }

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        PhysicsFunction<T_S,T_P>*                 solutionFunction,
        std::vector<Parameter*>                   parameters,
        const std::vector<unsigned int>&          order
        ) 
      : 
        m_parameters(parameters), 
        m_dimension( parameters.size() ),
        m_order(order)
    {
      m_solutionFunction = std::map< std::string, PhysicsFunction<T_S,T_P>* > ( 
          std::map<std::string, PhysicsFunction<T_S,T_P>* >(
            solutionFunction->name(),solutionFunction) );

      /* m_coefficients.resize( 1 ); */

      if (m_order.size() != parameters.size() )
      {
        std::cout 
          << std::endl
          << "\tERROR:"
          << " order vector dimension does not match number of parameters"
          << std::endl
          << std::endl;
        assert(0);
      }
    }

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        PhysicsFunction<T_S,T_P>*                 solutionFunction,
        std::vector<Parameter*>                   parameters,
        unsigned int                              order
        ) 
      : 
        m_parameters(parameters), 
        m_dimension( parameters.size() )
    {
      m_solutionFunction = std::map< std::string, PhysicsFunction<T_S,T_P>* > ( 
          std::map<std::string, PhysicsFunction<T_S,T_P>* >(
            solutionFunction->name(),solutionFunction) );

      m_order = std::vector<unsigned int>(m_dimension,order);
      /* m_coefficients.resize( 1 ); */
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

/********************************************//**
 * \brief Get the current expansion order
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<unsigned int> SurrogateModel<T_S,T_P>::getExpansionOrder( )
    const
  {
    return m_order;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P> 
    const std::map< std::string, std::vector<T_P> >
    SurrogateModel<T_S,T_P>::getCoefficients( ) const
    {
      return m_coefficients;
    }



}
#endif //SURROGATE_MODEL_H
