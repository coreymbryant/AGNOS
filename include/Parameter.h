
#ifndef PARAMETER_H
#define PARAMETER_H

#include <assert.h>
#include <boost/math/special_functions/legendre.hpp>

namespace AGNOS
{
  enum ParameterType { CONSTANT=0 , UNIFORM=1 };

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
      Parameter( int type, double min, double max );
      Parameter( std::string s, double min, double max );
      Parameter( );
      virtual ~Parameter( );

      // Manipulators
      int type() const; /**< 0:CONSTANT, 1:UNIFORM */
      double min() const;
      double max() const;
      double evalBasisPoly( int l, double x);


    private:

      int    m_type;
      double m_min;
      double m_max;


  };

/********************************************//**
 * \brief Simple constructor
 ***********************************************/
  Parameter::Parameter( int type, double min, double max )
    : m_type(type), m_min(min), m_max(max)
  { 
  }

/********************************************//**
 * \brief Simple constructor
 ***********************************************/
  Parameter::Parameter( std::string s, double min, double max )
    : m_min(min), m_max(max)
  { 
    if ( s == "CONSTANT" )
      m_type = CONSTANT;
    else if ( s == "UNIFORM" )
        m_type = UNIFORM;
    else
    {
      std::cout << "\n ERROR: Unrecognized parameter type\n\n" ;
      agnos_assert(0);
    }
  }

/********************************************//**
 * \brief 
 ***********************************************/
  Parameter::Parameter( )
    : m_type( 0 ), m_min(-1.0), m_max(1.0)
  { }

/********************************************//**
 * \brief 
 ***********************************************/
  Parameter::~Parameter( )
  {
  }

/********************************************//**
 * \brief 
 ***********************************************/
  int Parameter::type( ) const
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


    switch ( ParameterType(m_type) )
    {
      case CONSTANT:
        agnos_assert( (std::abs(m_min - m_max) <= 1e-16) );
        polyValue = 1.0 ;
        
        break;

      case UNIFORM:
        assert( (std::abs(m_min - m_max) > 1e-16) );

        scaledX = ( x - (m_min + m_max ) /2.0 ) * 2.0/(m_max-m_min) ;
        // 1/sqrt( 1/( 2*l +1 ) ) is normalization factor 1/sqrt(norm^2))
        polyValue = 1.0 / std::sqrt( 1.0/(2.0 * l + 1.0) )  
          * boost::math::legendre_p( l, scaledX);
        break;

      default:
        std::cout << "\n ERROR: Unrecognized parameter type\n\n" ;
        agnos_assert(0);

    }
    return polyValue;
  }

}

#endif // PARAMETER_H


