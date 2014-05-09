#ifndef PSEUDO_SPECTRAL_TENSOR_PRODUCT_H
#define PSEUDO_SPECTRAL_TENSOR_PRODUCT_H
#include "SurrogatePseudoSpectral.h"
#include "QuadratureTensorProduct.h"
#include <iostream>

namespace AGNOS
{

/********************************************//**
 * \brief Tensor product SurrogatePseudoSpectral class
 *
 * Tensor product quadrature is used to comptue coefficients in a
 * Pseudo-spectral method. Anisotropic polynomial order is supported
 * 
 ***********************************************/
  template<class T_S, class T_P>
    class PseudoSpectralTensorProduct : public SurrogatePseudoSpectral<T_S,T_P>
  {

    public:

      /** Constructor:  */
      PseudoSpectralTensorProduct( 
        const Communicator&               comm,
        std::shared_ptr<PhysicsModel<T_S,T_P> >               physics,
        const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
        const std::vector<unsigned int>&  order,
        std::set<std::string> computeSolutions = std::set<std::string>()
        );

      /** Secondary Constructor. 
       *  ***DO NOT MISTAKE FOR A COPY CONSTRUCTOR ***
       *  Intended use is for constructing a secondary surrogate model using the
       *  primary model as an evaluating object in the build routine. 
       *  If additional inputs are defined it will
       * construct a new surrogate increasing the order and using
       * primarySurrogate to perform evaluations in the constructions */
      PseudoSpectralTensorProduct( 
          std::shared_ptr<SurrogateModelBase<T_S,T_P> > primarySurrogate, 
          std::vector<unsigned int> increaseOrder = std::vector<unsigned int>(),
          unsigned int multiplyOrder = 1,
          std::set<std::string> evaluateSolutions = std::set<std::string>(),
          std::set<std::string> computeSolutions = std::set<std::string>()
          );

      /** Default destructor */
      ~PseudoSpectralTensorProduct( );

      /** Initialization routine */
      void initialize( ) ;

      /** return reference to quadrature rule */
      const QuadratureRule* getQuadRule( ) const ;


    protected:

      /** construct index set using recursion */
      void recurIndexSet(
          const unsigned int dim, 
          std::vector< std::vector<unsigned int> >& currentSet);

      QuadratureRule* _quadRule;


  };

  
}
#endif // PSEUDO_SPECTRAL_TENSOR_PRODUCT_H


