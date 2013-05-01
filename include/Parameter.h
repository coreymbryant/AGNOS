
#ifndef PARAMETER_H
#define PARAMETER_H

namespace AGNOS
{
  enum parameterType { UNIFORM = 0 };

  // TODO Need template ?
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


}

#endif // PARAMETER_H


