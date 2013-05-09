
#ifndef PARAMETER_H
#define PARAMETER_H

#include <assert.h>
#include <boost/math/special_functions/legendre.hpp>

namespace AGNOS
{
  enum parameterType { UNIFORM = 0 };

  /********************************************//**
   * \brief Uncertain parameter class
   *
   * Encapsulates the definition of uncertain/random parameters includes
   * parameter domain and type, as well as the appropriate orthogonal polynomial
   * definiton. 
   * 
   ***********************************************/
  class Parameter
  {

    public: 

      // constructor/destructor
      Parameter( enum parameterType, double min, double max );

      Parameter( );
      ~Parameter( );

      // Manipulators
      enum parameterType type() const; /**< 0: Uniform */
      double min() const;
      double max() const;
      double evalBasisPoly( int l, double x);

    private:

      enum parameterType(m_type);
      double m_min;
      double m_max;


  };

/********************************************//**
 * \brief Simple constructor
 ***********************************************/
  Parameter::Parameter( enum parameterType myType, double min, double max )
    : m_type(parameterType(myType)), m_min(min), m_max(max)
  { }

/********************************************//**
 * \brief 
 ***********************************************/
  Parameter::Parameter( )
    : m_type( UNIFORM ), m_min(-1.0), m_max(1.0)
  { }

/********************************************//**
 * \brief 
 ***********************************************/
  Parameter::~Parameter( )
  { }

/********************************************//**
 * \brief 
 ***********************************************/
  enum parameterType Parameter::type( ) const
  {
    return m_type;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  double Parameter::min( ) const
  {
    return m_min;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  double Parameter::max( ) const
  {
    return m_max;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  double Parameter::evalBasisPoly( int l, double x)
  {
    double polyValue;
    double scaledX; 

    switch ( m_type )
    {
      case UNIFORM:
        scaledX = ( x - (m_min + m_max ) /2.0 ) * 2.0/(m_max-m_min) ;
        // 1/sqrt( 1/( 2*l +1 ) ) is normalization factor 1/sqrt(norm^2))
        polyValue = 1.0 / std::sqrt( 1.0/(2.0 * l + 1.0) )  
          * boost::math::legendre_p( l, scaledX);
        break;

      default:
        std::cout << "\n ERROR: Unrecognized parameter type\n\n" ;
        assert(0);

    }
    return polyValue;
  }

}

#endif // PARAMETER_H


