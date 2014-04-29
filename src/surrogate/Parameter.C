
#include "Parameter.h"
#include <assert.h>
#include <boost/math/special_functions/legendre.hpp>

namespace AGNOS
{
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

/********************************************//**
 * \brief returns measure of [min,max] according to distribution 
 ***********************************************/
  double Parameter::measure( double min, double max )
  {
    double measure;

    // calculate measure based on distribution
    switch ( ParameterType(m_type) )
    {
      case CONSTANT:
        measure = 1.0 ;
        break;

      case UNIFORM:
        measure = (max-min)/(m_max-m_min);
        break;

      default:
        std::cout << "\n ERROR: Unrecognized parameter type\n\n" ;
        agnos_assert(0);

    }

    return measure;
  }

} // end namespace AGNOS

