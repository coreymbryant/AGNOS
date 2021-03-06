
#include "QuadratureTensorProduct.h"

// TensorProduct namespace since we may eventually have sparse grid as well
namespace AGNOS
{
/********************************************//**
 * \brief Rountine to populate quadPoints and quadWeights
 * 
 ***********************************************/
  QuadratureTensorProduct::QuadratureTensorProduct( 
      const std::vector<std::shared_ptr<AGNOS::Parameter> >& parameters,
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
 * \brief 
 * 
 ***********************************************/
  QuadratureTensorProduct::~QuadratureTensorProduct( )
  {}

/********************************************//**
 * \brief recurs to produce quad rule for higher dimensions
 *
 ***********************************************/
  void QuadratureTensorProduct::recurQuad(
      const std::vector<std::shared_ptr<AGNOS::Parameter> >& parameters,
      const int dim, const std::vector<unsigned int>& order, 
      double currentWeights[], double* currentPoints[] )
  {
    double* oneDimQuadWeights = new double[order[dim-1]+1];
    double* oneDimQuadPoints  = new double[order[dim-1]+1];
    oneDimQuadRule( *parameters[dim-1] , order[dim-1]+1, 
        oneDimQuadPoints, oneDimQuadWeights );

    if (dim == 1) 
    {
      // TODO change to assign instead of loop ?
      for(unsigned int id=0; id < order[dim-1]+1; id++)
      {
        currentWeights[id] = oneDimQuadWeights[id] ;
        currentPoints[id][0] = oneDimQuadPoints[id] ;
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
            currentWeights[outer] * oneDimQuadWeights[inner];
          
          for(int dir=0; dir < dim-1; dir++)
          {
            currentPoints[(outer)*(order[dim-1]+1) + inner][dir]
              = currentPoints[outer][dir];
          }

          currentPoints[(outer)*(order[dim-1]+1) + inner][dim-1] 
            = oneDimQuadPoints[inner];
        }
      }
    }
    delete[] oneDimQuadPoints;
    delete[] oneDimQuadWeights;

    return ;
  }

/********************************************//**
 * \brief construct 1D quad rule associated with parameter distribution
 * 
 * Currently only unifrom is supported
 ***********************************************/
    void QuadratureTensorProduct::oneDimQuadRule(
        const Parameter& parameter, const unsigned int n, 
        double oneDimQuadPoints[], double oneDimQuadWeights[] )
    {
      int myType = parameter.type() ;
      double min = parameter.min();
      double max = parameter.max();
      double scale, midpoint, scaleWeight;

      switch ( myType )
      {
        case CONSTANT:
          agnos_assert( (std::abs(min - max) <= 1e-16) );
          agnos_assert( (n == 1) );
          oneDimQuadWeights[0] = 1.;
          oneDimQuadPoints[0] = min;
          break;

        case UNIFORM:
          agnos_assert( (std::abs(min - max) > 1e-16) );
          // this class assumes uniform distribution
          webbur::legendre_compute( 
              n, oneDimQuadPoints, oneDimQuadWeights);
          
          // scale - scales to size of parameter domain
          scale = ( max - min)/2.0;
          // midpoint - shifts to midpoint of parameter domain
          midpoint = ( max + min )/2.0;

          for(unsigned int i=0; i < n; i++)
          {
            // 1/2 is to account for uniform weighting
            // scale * 1/(max-min) = 1/2
            oneDimQuadWeights[i] = oneDimQuadWeights[i] * 1./2.;
            
            // pt = midpoint + pt * scale
            oneDimQuadPoints[i] *= scale  ; 
            oneDimQuadPoints[i] += midpoint;
          }

          break;

        default:
          std::cout << "\n ERROR: Unrecognized parameter type\n\n" ;
          agnos_assert(0);
      }
        
    }


}

