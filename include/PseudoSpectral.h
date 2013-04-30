// TODO need to change double** into some sort of array of T_P

#ifndef PSEUDO_SPECTRAL_H
#define PSEUDO_SPECTRAL_H
#include "SurrogateModel.h"

namespace AGNOS
{

/********************************************//**
 * \brief Pseudo-spectral projection surrogate model
 *
 * This class provides the framework for constructing surrogate models using
 * non-intrusive spectral projection methods. 
 *
 * As of now it is only Tensor product quadrature to comptue coefficients but
 * could be extended to other methods as well. Both isotropic and non-isotropic
 * polynomial orders are supported. 
 *
 * Only Uniform Random variables at the moment
 * 
 ***********************************************/
  template<class T_S, class T_P>
    class PseudoSpectral : public SurrogateModel<T_S,T_P>
  {

    public:

      PseudoSpectral( 
          T_P* (*physicsFunction)(T_S& parameterValue),
          std::vector<Parameter*> parameters,
          unsigned int order 
          );
      PseudoSpectral( 
          T_P* (*physicsFunction)(T_S& parameterValue),
          std::vector<Parameter*> parameters,
          std::vector<unsigned int>& order
          );

      ~PseudoSpectral( );

      // these should carry over from surrogateModel right?
      /* virtual void build( ) = 0; */ 
      /* virtual T_P& evaluate( */ 
      /*     ) = 0; */
      /* virtual void refine( ) = 0; */

      // Manipulators
      void setExpansionOrder( std::vector<unsigned int>& order);
      std::vector<unsigned int> getExpansionOrder( ) const;

    private:

      std::vector<unsigned int> m_order;  // anisotropic

  };


/********************************************//**
 * \brief Constructor for isotropic order
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( 
        T_P* (*physicsFunction)(T_S& parameterValue),
        std::vector<Parameter*> parameters,
        unsigned int order
        )
      : SurrogateModel(physicsFunction,parameters)
    {
      m_order = std::vector<unsigned int>(m_dimension,order);
    }


/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( 
        T_P* (*physicsFunction)(T_S& parameterValue),
        std::vector<Parameter*> parameters,
        std::vector<unsigned int>& order
        )
      : SurrogateModel(physicsFunction,parameters), 
      m_order(order)
    {
    }

/********************************************//**
 * \brief Set expansion order
 ***********************************************/
  template<class T_S, class T_P>
    void PseudoSpectral<T_S,T_P>::setExpansionOrder( 
        std::vector<unsigned int>& order
        )
  {
    m_order = order ;
    return ;
  }

/********************************************//**
 * \brief Get the current expansion order
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<unsigned int> PseudoSpectral<T_S,T_P>::getExpansionOrder( )
    const
  {
    return m_order;
  }

  

  
}
#endif // PSEUDO_SPECTRAL_H


