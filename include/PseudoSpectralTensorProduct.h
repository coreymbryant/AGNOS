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
        const std::vector<unsigned int>&  order
        );

      /** Secondary Constructor. 
       *  ***DO NOT MISTAKE FOR A COPY CONSTRUCTOR ***
       *  Intended use is for constructing a secondary surrogate model using the
       *  primary model as an evaluating object in the build routine. 
       *  If additional inputs are defined it will
       * construct a new surrogate increasing the order and using
       * primarySurrogate to perform evaluations in the constructions */
      PseudoSpectralTensorProduct( 
          const SurrogateModel<T_S,T_P>* primarySurrogate, 
          unsigned int increaseOrder = 0,
          unsigned int multiplyOrder = 1,
          std::set<std::string> evaluateSolutions = std::set<std::string>(),
          std::set<std::string> computeSolutions = std::set<std::string>()
          );

      /** Default destructor */
      ~PseudoSpectralTensorProduct( );

      /** Initialization routine */
      void initialize( ) ;

      /** Refinement routine */
      void refine( );

      /** return reference to quadrature rule */
      const QuadratureRule* getQuadRule( ) const ;


    protected:

      /** construct index set using recursion */
      void recurIndexSet(
          const unsigned int dim, 
          std::vector< std::vector<unsigned int> >& currentSet);

      QuadratureRule* _quadRule;


  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectralTensorProduct<T_S,T_P>::PseudoSpectralTensorProduct( 
        const Communicator&               comm,
        std::shared_ptr<PhysicsModel<T_S,T_P> >               physics,
        const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
        const std::vector<unsigned int>& order
        )
      :
        SurrogatePseudoSpectral<T_S,T_P>( comm, physics, parameters, order)
    {
      initialize();
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectralTensorProduct<T_S,T_P>::PseudoSpectralTensorProduct( 
        const SurrogateModel<T_S,T_P>* primarySurrogate, 
        unsigned int increaseOrder ,
        unsigned int multiplyOrder ,
        std::set<std::string> evaluateSolutions,
        std::set<std::string> computeSolutions
        )
      : SurrogatePseudoSpectral<T_S,T_P>(primarySurrogate, increaseOrder,
          multiplyOrder, evaluateSolutions, computeSolutions)
    {
      initialize();
    }




/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectralTensorProduct<T_S,T_P>::~PseudoSpectralTensorProduct()
    {
      delete _quadRule;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void PseudoSpectralTensorProduct<T_S,T_P>::initialize( )
    {
      _quadRule = 
        new QuadratureTensorProduct( this->_parameters, this->_order);

      this->_nIntegrationPoints = ( _quadRule->getNQuadPoints() );

      this->_integrationPoints.clear();
      this->_integrationWeights.clear();
      this->_integrationPoints.resize(   this->_nIntegrationPoints );
      this->_integrationWeights.resize(  this->_nIntegrationPoints );

      double* quadWeights = _quadRule->getQuadWeights() ;
      this->_integrationWeights.assign( quadWeights,
          quadWeights+this->_nIntegrationPoints  );

      double** quadPoints = _quadRule->getQuadPoints() ;
      for (unsigned int point=0; point < this->_nIntegrationPoints; point++)
      {
        this->_integrationPoints[point] = T_P(this->_dimension);
        for (unsigned int dir=0; dir < this->_dimension; dir++)
          this->_integrationPoints[point](dir) = quadPoints[point][dir] ;
      }

      this->_indexSet.clear();
      this->_indexSet.reserve( this->_nIntegrationPoints );

      recurIndexSet( this->_dimension, this->_indexSet );

      this->_totalNCoeff = this->_indexSet.size();
      this->_coefficients.clear();
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P> 
    void PseudoSpectralTensorProduct<T_S,T_P>::recurIndexSet( 
          const unsigned int dim, 
          std::vector< std::vector<unsigned int> >& currentSet  )
    {

      if (dim == 1) 
      {
        for(unsigned int id=0; id < this->_order[dim-1]+1; id++)
          currentSet.push_back( std::vector<unsigned int>(1,id) );
      }

      else
      {
        recurIndexSet(dim-1,currentSet);

        int prevSize = currentSet.size();
        currentSet.resize( prevSize * (this->_order[dim-1]+1 ) );
        // keep track of how many array elements are non-zero
        /* int prevSize = order[0]+1; */
        /* for(unsigned int i=0; i < dim-2; i++) */
        /*   prevSize *= order[i+1] + 1; */

        for(int out=prevSize-1; out >= 0; out--)
        {
          for(int in=this->_order[dim-1]; in >= 0 ; in--)
          {
            currentSet[(out)*(this->_order[dim-1]+1) + in].reserve(dim);
            currentSet[(out)*(this->_order[dim-1]+1) + in] = currentSet[out];
            currentSet[(out)*(this->_order[dim-1]+1) + in].push_back(in);
          }
        }
      }


    }

  
  //TODO add refinement option to increase order a given amount instead of just
  //+1 in all directions

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void PseudoSpectralTensorProduct<T_S,T_P>::refine( )
    {
      for(unsigned int i=0; i<this->_dimension; i++)
      {
        this->_order[i]++;
      }
      this->initialize();
      this->build();
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P> 
    const QuadratureRule*
    PseudoSpectralTensorProduct<T_S,T_P>::getQuadRule( )
    const 
    {
      return this->_quadRule ;
    }






  
}
#endif // PSEUDO_SPECTRAL_TENSOR_PRODUCT_H


