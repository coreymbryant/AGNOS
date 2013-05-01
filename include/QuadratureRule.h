
#ifndef QUADRATURE_RULE_H
#define QUADRATURE_RULE_H

#include "sandia_rules.hpp"
#include "Parameter.h"

namespace AGNOS
{

  class QuadratureRule
  {
    public:

      QuadratureRule( );
      ~QuadratureRule( );

      // Manipulators
      unsigned int getDimension( ) const;
      unsigned int getNQuadPoints( ) const;
      double** getQuadPoints( ) const;
      double* getQuadWeights( ) const;
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
  QuadratureRule::QuadratureRule( )
  {
    // Do nothing, this needs to be defined for derived classes
  }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  QuadratureRule::~QuadratureRule( )
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
  unsigned int QuadratureRule::getNQuadPoints( ) const
  {
    return m_nQuadPoints;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  unsigned int QuadratureRule::getDimension( ) const
  {
    return m_dimension;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  double** QuadratureRule::getQuadPoints( ) const
  {
    return m_quadPoints;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  double* QuadratureRule::getQuadWeights( ) const
  {
    return m_quadWeights;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  void QuadratureRule::printQuadPoints( ) const
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
  void QuadratureRule::printQuadWeights( ) const
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
