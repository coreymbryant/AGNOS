
#ifndef PARAMETER_H
#define PARAMETER_H

#include "agnosDefines.h"
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


}

#endif // PARAMETER_H


