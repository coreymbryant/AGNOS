
#ifndef QUADRATURE_RULE_H
#define QUADRATURE_RULE_H

#include "Parameter.h"
#include "sandia_rules.hpp"

namespace AGNOS
{

  /********************************************//**
   * \brief Base class for Quadrature
   ***********************************************/
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


}

#endif // QUAD_RULE_H
