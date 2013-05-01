#ifndef TENSOR_PRODUCT_PSEUDO_SPECTRAL_H
#define TENSOR_PRODUCT_PSEUDO_SPECTRAL_H
#include "PseudoSpectral.h"
#include "TensorProductQuadrature.h"
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
      typedef T_P* (PhysicsModel<T_S,T_P>::*PhysicsMemberFcn)(T_S&);

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

      ~TensorProductPseudoSpectral( );

      // Surrogate constructors
      void build( ) ;
      T_P evaluate( 
          T_S& parameterValues /**< parameter values to evaluate*/
          );
      void refine( );

      // Manipulators
      const QuadratureRule<T_S>* getQuadRule( ) const ;


    protected:

      std::vector<T_P> m_coefficients;
      std::vector<double> m_indexSet;

      QuadratureRule<T_S>* m_quadRule;


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
      m_quadRule = new TensorProductQuadrature<T_S>( parameters, this->m_order);

      m_coefficients.resize( m_quadRule->getNQuadPoints() );

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
      m_quadRule = new TensorProductQuadrature<T_S>( parameters, this->m_order);

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
      /* T_P solVec; */ 

      /* unsigned int nQuadPoints = this->m_quadRule->getNQuadPoints(); */
      /* std::vector<T_S> parameterValues */
      /*   = this->m_quadRule->getQuadPoints(); */

      /* for (unsigned int j=0; j < nQuadPoints; j++) */
      /* { */
      /*   this->m_solutionFunction.compute(paramValues[j],solVec); */

      /*   for(unsigned int coeff=0; coeff < nQuadPoints; coeff++) */
      /*   { */
      /*     std::cout << "solVec[0] = " << solVec[0] << std::endl; */
      /*     std::cout << "poly = " */ 
      /*       << (this->m_parameters[0]->evalBasisPoly(coeff,paramVec[0]) ) */ 
      /*       << std::endl; */
      /*     std::cout << "weight = " */ 
      /*       << (this->m_quadRule->getQuadWeights())[0] */
      /*       << std::endl; */

      /*     for(unsigned int solId; solId < numberOfSolutionVectors; solId++) */
      /*       m_coefficients[coeff][solId][0] +=  solVec[solId][0] */
      /*         * (this->m_parameters[0]->evalBasisPoly(coeff,paramVec[0][0]) ) */ 
      /*         * (this->m_quadRule->getQuadWeights())[0] ; */
      /*   } */
      /*   std::cout << "coeff[" << j << "][0] = " << m_coefficients[j][0][0] << std::endl; */
      /* } */
     } 


/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    T_P TensorProductPseudoSpectral<T_S,T_P>::evaluate( 
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
    const QuadratureRule<T_S>*
    TensorProductPseudoSpectral<T_S,T_P>::getQuadRule( )
    const 
    {
      return this->m_quadRule ;
    }

  
}
#endif // TENSOR_PRODUCT_PSEUDO_SPECTRAL_H


