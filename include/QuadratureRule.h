
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


}

#endif // QUAD_RULE_H
