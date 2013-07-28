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
  template<class T_S>
    class PseudoSpectralTensorProduct : public SurrogatePseudoSpectral<T_S>
  {

    public:

      PseudoSpectralTensorProduct( 
          const Communicator&               comm,
          PhysicsModel<T_S>*                physics,
          const std::vector<Parameter*>     parameters,
          const unsigned int                order 
          );
      PseudoSpectralTensorProduct( 
          const Communicator&               comm,
          PhysicsModel<T_S>*                physics,
          const std::vector<Parameter*>     parameters,
          const std::vector<unsigned int>&  order
          );

      void initialize( ) ;

      ~PseudoSpectralTensorProduct( );

      void refine( );

      // Manipulators
      const QuadratureRule* getQuadRule( ) const ;


    protected:

      void recurIndexSet(
          const unsigned int dim, 
          std::vector< std::vector<unsigned int> >& currentSet);

      QuadratureRule* _quadRule;


  };


/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S>
    PseudoSpectralTensorProduct<T_S>::PseudoSpectralTensorProduct( 
        const Communicator&               comm,
        PhysicsModel<T_S>*                physics,
        const std::vector<Parameter*> parameters,
        const unsigned int order
        )
      : SurrogatePseudoSpectral<T_S>(comm,physics,parameters,order)
    {
      initialize();
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S>
    PseudoSpectralTensorProduct<T_S>::PseudoSpectralTensorProduct( 
        const Communicator&               comm,
        PhysicsModel<T_S>*                physics,
        const std::vector<Parameter*> parameters,
        const std::vector<unsigned int>& order
        )
      : SurrogatePseudoSpectral<T_S>(comm,physics,parameters,order)
    {
      initialize();
    }


/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S>
    PseudoSpectralTensorProduct<T_S>::~PseudoSpectralTensorProduct()
    {
      delete _quadRule;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S>
    void PseudoSpectralTensorProduct<T_S>::initialize( )
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
        this->_integrationPoints[point] = std::vector<Number>(this->_dimension);
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
  template<class T_S> 
    void PseudoSpectralTensorProduct<T_S>::recurIndexSet( 
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
  template<class T_S>
    void PseudoSpectralTensorProduct<T_S>::refine( )
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
  template<class T_S> 
    const QuadratureRule*
    PseudoSpectralTensorProduct<T_S>::getQuadRule( )
    const 
    {
      return this->_quadRule ;
    }






  
}
#endif // PSEUDO_SPECTRAL_TENSOR_PRODUCT_H


