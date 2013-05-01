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

      typedef T_P* (PhysicsModel<T_S,T_P>::*PhysicsMemberFcn)(T_S&);
      PseudoSpectral( 
          PhysicsFunction<T_S,T_P>& solutionFunction,
          const std::vector<Parameter*> parameters,
          const unsigned int order 
          );
      PseudoSpectral( 
          PhysicsFunction<T_S,T_P>& solutionFunction,
          const std::vector<Parameter*> parameters,
          const std::vector<unsigned int>& order
          );

      virtual ~PseudoSpectral( );

      // these should carry over from surrogateModel right?
      /* virtual void build( ) = 0; */ 
      /* virtual T_P& evaluate( */ 
      /*     ) = 0; */
      /* virtual void refine( ) = 0; */

      // Manipulators
      void setExpansionOrder( std::vector<unsigned int>& order);
      std::vector<unsigned int> getExpansionOrder( ) const;

    protected:

      std::vector<unsigned int> m_order;  // anisotropic

  };


/********************************************//**
 * \brief Constructor for isotropic order
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( 
        PhysicsFunction<T_S,T_P>& solutionFunction,
        const std::vector<Parameter*> parameters,
        const unsigned int order
        )
      : SurrogateModel<T_S,T_P>(solutionFunction,parameters)
    {
      m_order = std::vector<unsigned int>(this->m_dimension,order);
    }


/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( 
        PhysicsFunction<T_S,T_P>& solutionFunction,
        const std::vector<Parameter*> parameters,
        const std::vector<unsigned int>& order
        )
      : SurrogateModel<T_S,T_P>(solutionFunction,parameters), 
      m_order(order)
    {
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

/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::~PseudoSpectral( )
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


