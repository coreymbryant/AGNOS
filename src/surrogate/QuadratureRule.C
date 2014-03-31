
#include "QuadratureRule.h"

namespace AGNOS
{
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

