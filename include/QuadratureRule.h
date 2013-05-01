
#ifndef QUADRATURE_RULE_H
#define QUADRATURE_RULE_H

#include "sandia_rules.hpp"
#include "Parameter.h"

namespace AGNOS
{

  template<class T_S>
  class QuadratureRule
  {
    public:

      QuadratureRule( );
      ~QuadratureRule( );

      // Manipulators
      unsigned int getDimension( ) const;
      unsigned int getNQuadPoints( ) const;
      std::vector<T_S> getQuadPoints( ) const;
      T_S getQuadWeights( ) const;
      void printQuadWeights( ) const;
      void printQuadPoints( ) const;

    protected:
      unsigned int m_dimension;
      unsigned int m_nQuadPoints;
      double** m_quadPoints;
      double* m_quadWeights;
  };

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  template<class T_S>
  QuadratureRule<T_S>::QuadratureRule( )
  {
    // Do nothing, this needs to be defined for derived classes
  }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  template<class T_S>
  QuadratureRule<T_S>::~QuadratureRule( )
  {
    for (unsigned int i=0; i< m_nQuadPoints; i++)
      delete [] m_quadPoints[i] ;

    delete [] m_quadPoints;
    delete [] m_quadWeights;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S>
  unsigned int QuadratureRule<T_S>::getNQuadPoints( ) const
  {
    return m_nQuadPoints;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S>
  unsigned int QuadratureRule<T_S>::getDimension( ) const
  {
    return m_dimension;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S>
  std::vector<T_S> QuadratureRule<T_S>::getQuadPoints( ) const
  {
    std::vector<T_S> quadPoints;
    quadPoints.resize(m_nQuadPoints);
    for (unsigned int n=0; n< m_nQuadPoints; n++)
    {
      quadPoints[n] = T_S(m_dimension,0.0);
      for(unsigned int d=0; d< m_dimension; d++)
        quadPoints[n][d] = m_quadPoints[n][d];
    }
    return quadPoints;
  }
/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S>
  T_S QuadratureRule<T_S>::getQuadWeights( ) const
  {
    std::vector<T_S> quadWeights;
    quadWeights.resize(m_nQuadPoints);
    for (unsigned int n=0; n< m_nQuadPoints; n++)
      quadWeights[n] = m_quadWeights[n];

    return quadWeights;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S>
  void QuadratureRule<T_S>::printQuadPoints( ) const
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
  template<class T_S>
  void QuadratureRule<T_S>::printQuadWeights( ) const
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

}

#endif // QUAD_RULE_H
