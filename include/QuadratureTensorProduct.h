
#ifndef QUADRATURE_TENSOR_PRODUCT_H
#define QUADRATURE_TENSOR_PRODUCT_H
#include <assert.h>
#include "QuadratureRule.h"

// TensorProduct namespace since we may eventually have sparse grid as well
namespace AGNOS
{
  /********************************************//**
   * \brief Tensor product QuadratureRule
   *
   * Constructs QuadratureRule for a tensor product of Parameters.
   * Supports anisotropy
   *
   * 
   ***********************************************/
  class QuadratureTensorProduct : public QuadratureRule
  {

    public:

      QuadratureTensorProduct( 
          const std::vector<std::shared_ptr<AGNOS::Parameter> >& parameters,
          const std::vector<unsigned int>& order
          );

      QuadratureTensorProduct( );
      ~QuadratureTensorProduct( );

    protected:

      void recurQuad(
          const std::vector<std::shared_ptr<Parameter> >& parameters, 
          const int dim, const std::vector<unsigned int>& order, 
          double currentWeights[], double* currentPoints[] );

      void oneDimQuadRule(
        const Parameter& parameter, const unsigned int n, 
        double oneDimQuadPoints[], double oneDimQuadWeights[] );

  };


}

#endif // QUADRATURE_TENSOR_PRODUCT_H

