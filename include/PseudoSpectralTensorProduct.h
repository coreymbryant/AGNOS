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

      PseudoSpectralTensorProduct( 
          const Communicator*               comm,
          PhysicsFunction<T_S,T_P>*         solutionFunction,
          const std::vector<Parameter*>     parameters,
          const unsigned int                order 
          );
      PseudoSpectralTensorProduct( 
          const Communicator*               comm,
          PhysicsFunction<T_S,T_P>*         solutionFunction,
          const std::vector<Parameter*>     parameters,
          const std::vector<unsigned int>&  order
          );

      PseudoSpectralTensorProduct( 
          const Communicator*               comm,
          std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
          const std::vector<Parameter*>                       parameters,
          const unsigned int                                  order 
          );
      PseudoSpectralTensorProduct( 
          const Communicator*               comm,
          std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
          const std::vector<Parameter*>                       parameters,
          const std::vector<unsigned int>&                    order
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

      QuadratureRule* m_quadRule;


  };


/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectralTensorProduct<T_S,T_P>::PseudoSpectralTensorProduct( 
        const Communicator*               comm,
        PhysicsFunction<T_S,T_P>* solutionFunction,
        const std::vector<Parameter*> parameters,
        const unsigned int order
        )
      : SurrogatePseudoSpectral<T_S,T_P>(comm,solutionFunction,parameters,order)
    {
      initialize();
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectralTensorProduct<T_S,T_P>::PseudoSpectralTensorProduct( 
        const Communicator*               comm,
        PhysicsFunction<T_S,T_P>* solutionFunction,
        const std::vector<Parameter*> parameters,
        const std::vector<unsigned int>& order
        )
      : SurrogatePseudoSpectral<T_S,T_P>(comm,solutionFunction,parameters,order)
    {
      initialize();
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectralTensorProduct<T_S,T_P>::PseudoSpectralTensorProduct( 
        const Communicator*               comm,
        std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
        const std::vector<Parameter*>                       parameters,
        const unsigned int                                  order
        )
      : SurrogatePseudoSpectral<T_S,T_P>(comm,solutionFunction,parameters,order)
    {
      initialize();
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectralTensorProduct<T_S,T_P>::PseudoSpectralTensorProduct( 
        const Communicator*               comm,
        std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
        const std::vector<Parameter*>                       parameters,
        const std::vector<unsigned int>&                    order
        )
      : SurrogatePseudoSpectral<T_S,T_P>(comm,solutionFunction,parameters,order)
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
      delete m_quadRule;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void PseudoSpectralTensorProduct<T_S,T_P>::initialize( )
    {
      m_quadRule = 
        new QuadratureTensorProduct( this->m_parameters, this->m_order);

      this->m_nIntegrationPoints = ( m_quadRule->getNQuadPoints() );

      this->m_integrationPoints.clear();
      this->m_integrationWeights.clear();
      this->m_integrationPoints.resize(   this->m_nIntegrationPoints );
      this->m_integrationWeights.resize(  this->m_nIntegrationPoints );

      double* quadWeights = m_quadRule->getQuadWeights() ;
      this->m_integrationWeights.assign( quadWeights,
          quadWeights+this->m_nIntegrationPoints  );

      double** quadPoints = m_quadRule->getQuadPoints() ;
      for (unsigned int point=0; point < this->m_nIntegrationPoints; point++)
      {
        this->m_integrationPoints[point] = T_P(this->m_dimension);
        for (unsigned int dir=0; dir < this->m_dimension; dir++)
          this->m_integrationPoints[point](dir) = quadPoints[point][dir] ;
      }

      this->m_indexSet.clear();
      this->m_indexSet.reserve( this->m_nIntegrationPoints );

      recurIndexSet( this->m_dimension, this->m_indexSet );

      this->m_totalNCoeff = this->m_indexSet.size();
      this->m_coefficients.clear();
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
        for(unsigned int id=0; id < this->m_order[dim-1]+1; id++)
          currentSet.push_back( std::vector<unsigned int>(1,id) );
      }

      else
      {
        recurIndexSet(dim-1,currentSet);

        int prevSize = currentSet.size();
        currentSet.resize( prevSize * (this->m_order[dim-1]+1 ) );
        // keep track of how many array elements are non-zero
        /* int prevSize = order[0]+1; */
        /* for(unsigned int i=0; i < dim-2; i++) */
        /*   prevSize *= order[i+1] + 1; */

        for(int out=prevSize-1; out >= 0; out--)
        {
          for(int in=this->m_order[dim-1]; in >= 0 ; in--)
          {
            currentSet[(out)*(this->m_order[dim-1]+1) + in].reserve(dim);
            currentSet[(out)*(this->m_order[dim-1]+1) + in] = currentSet[out];
            currentSet[(out)*(this->m_order[dim-1]+1) + in].push_back(in);
          }
        }
      }


    }

  

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void PseudoSpectralTensorProduct<T_S,T_P>::refine( )
    {
      for(unsigned int i=0; i<this->m_dimension; i++)
      {
        this->m_order[i]++;
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
      return this->m_quadRule ;
    }






  
}
#endif // PSEUDO_SPECTRAL_TENSOR_PRODUCT_H


