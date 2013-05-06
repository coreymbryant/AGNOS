
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
          const std::vector<Parameter*>& parameters,
          const std::vector<unsigned int>& order
          );

      QuadratureTensorProduct( );
      ~QuadratureTensorProduct( );

    protected:

      void recurQuad(
          const std::vector<Parameter*>& parameters,
          const int dim, const std::vector<unsigned int>& order, 
          double currentWeights[], double* currentPoints[] );

      void oneDimQuadRule(
        enum parameterType, const unsigned int order, 
        double oneDimQuadPoints[], double oneDimQuadWeights[] );

  };


/********************************************//**
 * \brief Rountine to populate quadPoints and quadWeights
 * 
 ***********************************************/
  QuadratureTensorProduct::QuadratureTensorProduct( 
      const std::vector<Parameter*>& parameters, 
      const std::vector<unsigned int>& order 
      )
  {
    this->m_dimension = parameters.size();

    this->m_nQuadPoints = order[0]+1;
    for(unsigned int i=1; i< this->m_dimension; i++)
      this->m_nQuadPoints *= (order[i]+1);

    this->m_quadPoints = new double*[this->m_nQuadPoints];
    for(unsigned int i=0; i< this->m_nQuadPoints; i++)
      this->m_quadPoints[i] = new double[this->m_dimension];
    this->m_quadWeights = new double[this->m_nQuadPoints];

    recurQuad(
        parameters,this->m_dimension,order,
        this->m_quadWeights,this->m_quadPoints);

    return ;
  }

/********************************************//**
 * \brief recurs to produce quad rule for higher dimensions
 *
 ***********************************************/
  void QuadratureTensorProduct::recurQuad(
      const std::vector<Parameter*>& parameters,
      const int dim, const std::vector<unsigned int>& order, 
      double currentWeights[], double* currentPoints[] )
  {
    double* oneDimQuadWeights = new double[order[dim-1]+1];
    double* oneDimQuadPoints = new double[order[dim-1]+1];
    oneDimQuadRule( (parameters[dim-1])->type() , order[dim-1]+1, 
        oneDimQuadPoints, oneDimQuadWeights );

    double scaling = 
      ( (parameters[dim-1])->max() - (parameters[dim-1])->min())/2.0;
    double midpoint = 
      ( (parameters[dim-1])->max() + (parameters[dim-1])->min() )/2.0;

    if (dim == 1) 
    {
      for(unsigned int id=0; id < order[dim-1]+1; id++)
      {
        /* currentPoints[dim-1][id] = midpoint */ 
        /*   + oneDimQuadPoints[id] * scaling; */
          
        currentWeights[id] = oneDimQuadWeights[id] * scaling;
        currentPoints[id][0] = midpoint 
          + oneDimQuadPoints[id] * scaling;
      }
    }

    else 
    {
      recurQuad(parameters,dim-1,order,currentWeights,currentPoints);
      
      // keep track of how many array elements are non-zero
      int prevSize = order[0]+1;
      for(unsigned int i=0; i < dim-2; i++)
        prevSize *= order[i+1] + 1;

      for(int outer=prevSize-1; outer >= 0; outer--)
      {
        for(int inner=order[dim-1]; inner >= 0 ; inner--)
        {
          currentWeights[(outer)*(order[dim-1]+1) + inner] = 
            currentWeights[outer] * oneDimQuadWeights[inner] * scaling;
          
          for(int dir=0; dir < dim-1; dir++)
          {
            currentPoints[(outer)*(order[dim-1]+1) + inner][dir]
              = currentPoints[outer][dir];
          }

          currentPoints[(outer)*(order[dim-1]+1) + inner][dim-1] 
            = midpoint + oneDimQuadPoints[inner] * scaling;
        }
      }
    }
    return ;
  }

/********************************************//**
 * \brief construct 1D quad rule associated with parameter distribution
 * 
 * Currently only unifrom is supported
 ***********************************************/
    void QuadratureTensorProduct::oneDimQuadRule(
        enum parameterType myType, const unsigned int order, 
        double oneDimQuadPoints[], double oneDimQuadWeights[] )
    {
      switch ( myType )
      {
        case UNIFORM:
          // this class assumes uniform distribution
          webbur::legendre_compute( 
              order, oneDimQuadPoints, oneDimQuadWeights);
          break;

        default:
          std::cout << "\n ERROR: Unrecognized parameter type\n\n" ;
          assert(0);
      }
        
    }


}

#endif // QUADRATURE_TENSOR_PRODUCT_H

