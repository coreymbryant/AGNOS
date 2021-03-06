
#include "PseudoSpectralTensorProduct.h"
#include <iostream>

namespace AGNOS
{
/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectralTensorProduct<T_S,T_P>::PseudoSpectralTensorProduct( 
        const Communicator&               comm,
        std::shared_ptr<PhysicsModel<T_S,T_P> >               physics,
        const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
        const std::vector<unsigned int>& order,
        std::set<std::string> computeSolutions
        )
      :
        SurrogateModelBase<T_S,T_P>(comm,parameters,order,computeSolutions),
        SurrogatePseudoSpectral<T_S,T_P>( comm, physics,
            parameters, order, computeSolutions)
    {
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectralTensorProduct<T_S,T_P>::PseudoSpectralTensorProduct( 
        std::shared_ptr<SurrogateModelBase<T_S,T_P> > primarySurrogate, 
        std::vector<unsigned int> increaseOrder ,
        unsigned int multiplyOrder ,
        std::set<std::string> evaluateSolutions,
        std::set<std::string> computeSolutions
        )
        : 
          SurrogateModelBase<T_S,T_P>(
            primarySurrogate->getComm(),
            primarySurrogate->getParameters(),
            primarySurrogate->getExpansionOrder(),
            computeSolutions),
         SurrogatePseudoSpectral<T_S,T_P>(primarySurrogate, increaseOrder,
            multiplyOrder, evaluateSolutions, computeSolutions)
    {
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
    const QuadratureRule*
    PseudoSpectralTensorProduct<T_S,T_P>::getQuadRule( )
    const 
    {
      return this->_quadRule ;
    }



  template class PseudoSpectralTensorProduct<libMesh::DenseVector<double>,
           libMesh::DenseVector<double> >;




}
