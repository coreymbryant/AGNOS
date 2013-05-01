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
          std::vector<Parameter*> parameters,
          unsigned int order 
          );
      TensorProductPseudoSpectral( 
          PhysicsMemberFcn physicsFunction,
          std::vector<Parameter*> parameters,
          std::vector<unsigned int>& order
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
      void printQuadWeights( ) const;
      void printQuadPoints( ) const;

    protected:

      std::vector<T_P> m_coefficients;
      std::vector<double> m_indexSet;

      unsigned int m_nQuadPoints;
      double**  m_quadPoints;
      double* m_quadWeights;

      void constructQuadRule();
      

  };


/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    TensorProductPseudoSpectral<T_S,T_P>::TensorProductPseudoSpectral( 
        PhysicsMemberFcn physicsFunction,
        std::vector<Parameter*> parameters,
        unsigned int order
        )
      : PseudoSpectral<T_S,T_P>(physicsFunction,parameters,order)
    {
      m_nQuadPoints = this->m_order[0]+1;
      for(unsigned int i=1; i< this->m_dimension; i++)
        m_nQuadPoints *= (this->m_order[i]+1);

      m_coefficients.resize(m_nQuadPoints);

      m_quadPoints = new double*[m_nQuadPoints];
      for(unsigned int i=0; i<m_nQuadPoints; i++)
        m_quadPoints[i] = new double[this->m_dimension];
      m_quadWeights = new double[m_nQuadPoints];

      constructQuadRule();
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    TensorProductPseudoSpectral<T_S,T_P>::TensorProductPseudoSpectral( 
        PhysicsMemberFcn physicsFunction,
        std::vector<Parameter*> parameters,
        std::vector<unsigned int>& order
        )
      : PseudoSpectral<T_S,T_P>(physicsFunction,parameters,order)
    {
      m_nQuadPoints = this->m_order[0]+1;
      for(unsigned int i=1; i< this->m_dimension; i++)
        m_nQuadPoints *= (this->m_order[i]+1);

      m_coefficients.resize(m_nQuadPoints);

      m_quadPoints = new double*[m_nQuadPoints];
      for(unsigned int i=0; i<m_nQuadPoints; i++)
        m_quadPoints[i] = new double[this->m_dimension];
      m_quadWeights = new double[m_nQuadPoints];

      constructQuadRule();
    }


  
/********************************************//**
 * \brief Rountine to populate quadPoints and quadWeights
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
    void TensorProductPseudoSpectral<T_S,T_P>::constructQuadRule( )
    {

      TensorProduct::constructQuadRule(
          this->m_parameters,
          this->m_order,
          m_quadWeights,
          m_quadPoints
          );

      return ;
    }
  

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    TensorProductPseudoSpectral<T_S,T_P>::~TensorProductPseudoSpectral()
    {
      for (unsigned int i=0; i < m_nQuadPoints; i++)
        delete [] m_quadPoints[i];
      delete [] m_quadWeights;
    }


/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void TensorProductPseudoSpectral<T_S,T_P>::printQuadPoints( ) const
    {
      std::cout << std::endl;
      std::cout << "====================================================" <<
        std::endl;
      std::cout << " Quadrature points " << std::endl;
      std::cout << "----------------------------------------------------" <<
        std::endl;
      std::cout << "  \\ x        " << std::endl;
      std::cout << "   \\" ;
      for(unsigned int dim=0; dim < this->m_dimension; dim++)
        std::cout << std::setw(12) << "x_" << dim << " " ;
      std::cout << std::endl;
      std::cout << " id \\  " << std::endl;
      std::cout << "----------------------------------------------------" <<
        std::endl;
      for(int ix=0; ix < m_nQuadPoints ; ix++)  
      {
        std::cout << std::setw(3) << ix << "  |  " ;
        for(int iy=0; iy < this->m_dimension; iy++)
        {
          std::cout << std::scientific << std::setprecision(5) << std::setw(12)
            << m_quadPoints[ix][iy] << "  ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
      return;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
    void TensorProductPseudoSpectral<T_S,T_P>::printQuadWeights( ) const
    {
      double sum = 0.0;
      std::cout << std::endl;
      std::cout << "====================================================" <<
        std::endl;
      std::cout << " Quadrature weights " << std::endl;
      std::cout << "----------------------------------------------------" <<
        std::endl;
      for(int ip=0; ip < m_nQuadPoints; ip++){  
        std::cout << ip << "   | ";
        std::cout << m_quadWeights[ip] << "  ";
        std::cout << std::endl;
        sum += m_quadWeights[ip];
      }
      std::cout << "Sum = " << sum << std::endl;
      std::cout << std::endl;
      return;
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
      for(unsigned int pt=0; pt < m_nQuadPoints; pt++)
      {
        T_P ptSol = this->m_physicsFunction(m_quadPoints[pt]);
        for(unsigned int coeff; coeff < m_nQuadPoints; coeff++)
          // coeff += sol * poly * weight
          m_coefficients[coeff] += ptSol;
          /*   * POLY[coeff] * m_quadWeights[coeff] ; */
      }
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


  
}
#endif // TENSOR_PRODUCT_PSEUDO_SPECTRAL_H


