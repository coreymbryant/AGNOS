#ifndef TENSOR_PRODUCT_PSEUDO_SPECTRAL_H
#define TENSOR_PRODUCT_PSEUDO_SPECTRAL_H
#include "SurrogateModelPseudoSpectral.h"
#include "QuadratureRuleTensorProduct.h"
#include <iostream>

namespace AGNOS
{

/********************************************//**
 * \brief Tensor product Pseudo-spectral projection 
 *
 * Tensor product quadrature is used to comptue coefficients in a
 * Pseudo-spectral method
 * 
 ***********************************************/
  template<class T_S, class T_P>
    class TensorProductPseudoSpectral : public PseudoSpectral<T_S,T_P>
  {

    public:

      TensorProductPseudoSpectral( 
        PhysicsFunction<T_S,T_P>& solutionFunction,
          const std::vector<Parameter*> parameters,
          const unsigned int order 
          );
      TensorProductPseudoSpectral( 
        PhysicsFunction<T_S,T_P>& solutionFunction,
          const std::vector<Parameter*> parameters,
          const std::vector<unsigned int>& order
          );
      void initialize( ) ;
      void recurIndexSet( 
          );

      ~TensorProductPseudoSpectral( );

      void refine( );

      // Manipulators
      const QuadratureRule* getQuadRule( ) const ;


    protected:

      QuadratureRule* m_quadRule;


  };


/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    TensorProductPseudoSpectral<T_S,T_P>::TensorProductPseudoSpectral( 
        PhysicsFunction<T_S,T_P>& solutionFunction,
        const std::vector<Parameter*> parameters,
        const unsigned int order
        )
      : PseudoSpectral<T_S,T_P>(solutionFunction,parameters,order)
    {
      initialize();
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    TensorProductPseudoSpectral<T_S,T_P>::TensorProductPseudoSpectral( 
        PhysicsFunction<T_S,T_P>& solutionFunction,
        const std::vector<Parameter*> parameters,
        const std::vector<unsigned int>& order
        )
      : PseudoSpectral<T_S,T_P>(solutionFunction,parameters,order)
    {
      initialize();
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void TensorProductPseudoSpectral<T_S,T_P>::initialize( )
    {
      m_quadRule = 
        new TensorProductQuadrature( this->m_parameters, this->m_order);

      this->m_nIntegrationPoints = ( m_quadRule->getNQuadPoints() );

      this->m_integrationPoints.resize(   this->m_nIntegrationPoints );
      this->m_integrationWeights.resize(  this->m_nIntegrationPoints );

      double* quadWeights = m_quadRule->getQuadWeights() ;
      this->m_integrationWeights.assign( quadWeights,
          quadWeights+this->m_nIntegrationPoints  );

      double** quadPoints = m_quadRule->getQuadPoints() ;
      for (unsigned int point=0; point < this->m_nIntegrationPoints; point++)
      {
        this->m_integrationPoints[point] = T_P(this->m_dimension);
        for (unsigned int dir=0; dir < this->m_dimension; dir++)
          this->m_integrationPoints[point](dir) = quadPoints[point][dir] ;
      }

      this->m_coefficients.resize( this->m_nIntegrationPoints );
      this->m_indexSet.resize( this->m_nIntegrationPoints );
      for (unsigned int id=0; id < this->m_indexSet.size(); id++)
      {
        this->m_indexSet[id].resize( this->m_dimension );
        for (unsigned int dir=0; dir< this->m_dimension; dir++)
          for (unsigned int ord=0; ord < this->m_order[dir]; ord++)
            this->m_indexSet[id][dir] = 

      }
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P> 
    TensorProductPseudoSpectral<T_S,T_P>::recurIndexSet( 
        )
    {
    }

  

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    TensorProductPseudoSpectral<T_S,T_P>::~TensorProductPseudoSpectral()
    {
      delete m_quadRule;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P> 
    const QuadratureRule*
    TensorProductPseudoSpectral<T_S,T_P>::getQuadRule( )
    const 
    {
      return this->m_quadRule ;
    }





/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void TensorProductPseudoSpectral<T_S,T_P>::refine( )
    {
      /* for(unsigned int i=0; i<m_dimension; i++) */
      /*   m_order[i]++; */
    }

  
}
#endif // TENSOR_PRODUCT_PSEUDO_SPECTRAL_H


