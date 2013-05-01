#ifndef TENSOR_PRODUCT_PSEUDO_SPECTRAL_H
#define TENSOR_PRODUCT_PSEUDO_SPECTRAL_H
#include "PseudoSpectral.h"
#include "TensorProductQuadrature.h"

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
      typedef T_P* (PhysicsModel<T_S,T_P>::*PhysicsMemberFcn)(T_S&);

      TensorProductPseudoSpectral( 
          PhysicsMemberFcn physicsFunction,
          const std::vector<Parameter*> parameters,
          const unsigned int order 
          );
      TensorProductPseudoSpectral( 
          PhysicsMemberFcn physicsFunction,
          const std::vector<Parameter*> parameters,
          const std::vector<unsigned int>& order
          );

      ~TensorProductPseudoSpectral( );

      // Surrogate constructors
      // TODO
      void build( ) ;
      T_P* evaluate( 
          T_S& parameterValues /**< parameter values to evaluate*/
          );
      void refine( );

      // Manipulators
      const QuadratureRule* getQuadRule( ) const ;


    protected:

      std::vector<T_P> m_coefficients;
      std::vector<double> m_indexSet;

      QuadratureRule* m_quadRule;


  };


/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    TensorProductPseudoSpectral<T_S,T_P>::TensorProductPseudoSpectral( 
        PhysicsMemberFcn physicsFunction,
        const std::vector<Parameter*> parameters,
        const unsigned int order
        )
      : PseudoSpectral<T_S,T_P>(physicsFunction,parameters,order)
    {
      m_quadRule = new TensorProductQuadrature( parameters, this->m_order);
      m_coefficients.resize( m_quadRule->getNQuadPoints() );
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    TensorProductPseudoSpectral<T_S,T_P>::TensorProductPseudoSpectral( 
        PhysicsMemberFcn physicsFunction,
        const std::vector<Parameter*> parameters,
        const std::vector<unsigned int>& order
        )
      : PseudoSpectral<T_S,T_P>(physicsFunction,parameters,order)
    {
      m_quadRule = new TensorProductQuadrature( parameters, this->m_order);
      m_coefficients.resize( m_quadRule->getNQuadPoints() );
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
    void TensorProductPseudoSpectral<T_S,T_P>::build( )
    {
      // TODO
      /* for(unsigned int pt=0; pt < m_nQuadPoints; pt++) */
      /* { */
      /*   T_P ptSol = this->m_physicsFunction(m_quadPoints[pt]); */
      /*   for(unsigned int coeff; coeff < m_nQuadPoints; coeff++) */
      /*     // coeff += sol * poly * weight */
      /*     m_coefficients[coeff] += ptSol; */
      /*     /1*   * POLY[coeff] * m_quadWeights[coeff] ; *1/ */
      /* } */
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    T_P* TensorProductPseudoSpectral<T_S,T_P>::evaluate( 
        T_S& parameterValues /**< parameter values to evaluate*/
        )
    {
      /* T_P currSol; */
      /* // TODO */
      /* for(unsigned int coeff=0; coeff < m_nQuadPoints; coeff++) */
      /*   /1* currSol += coeff * Poly ; *1/ */
      /*   currSol = 0; */

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

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P> 
    const QuadratureRule* TensorProductPseudoSpectral<T_S,T_P>::getQuadRule( )
    const 
    {
      return this->m_quadRule ;
    }

  
}
#endif // TENSOR_PRODUCT_PSEUDO_SPECTRAL_H


