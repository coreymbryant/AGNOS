
#ifndef PSEUDO_SPECTRAL_H
#define PSEUDO_SPECTRAL_H

#include "SurrogateModel.h"
#include "sandia_rules.hpp"

namespace AGNOS
{

  template<ParameterDataType, PhysicsDataType>
    class PseudoSpectral : public SurrogateModel
  {

    public:

      PseudoSpectral( );
      PseudoSpectral( 
          unsigned int order, 
          unsigned int dimension, 
          std::vector<double>& m_mins, 
          std::vector<double>& m_maxs
          );
      ~PseudoSpectral( );

    private:

      unsigned int m_order;       // isotropic for now
      unsigned int m_dimension;   // all uniform distributed
      std::vector<double> m_mins; // can have different range though
      std::vector<double> m_maxs;

      std::vector<double> m_coefficients;
      std::vecotr<double> m_indexSet;

      std::vector< std::vector<double> >  m_quadPoints;
      std::vector<double> m_quadWeights;

  };


  template<ParameterDataType, PhysicsDataType>
    PseudoSpectral::PseudoSpectral( )
      : m_order = 0, dimension = 1
    {
      m_mins.push_back(-1.0)
      m_maxs.push_back(1.0);

      points = new double[1];
      weights = new double[1];

      webbur::legendre_compute( 1, points, weights );

      m_quadPoints.push_back( points );
      m_quadWeights = *weights;

      return;
    }

  template<ParameterDataType, PhysicsDataType>
    PseudoSpectral::PseudoSpectral( 
        unsigned int order, 
        unsigned int dimension, 
        std::vector<double>& m_mins, 
        std::vector<double>& m_maxs
        )
    {
      PseudoSpectral( );
      return;
    }

  template<ParameterDataType, PhysicsDataType>
    PseudoSpectral::~PseudoSpectral()
    {
      return ;
    }


  template<ParameterDataType, PhysicsDataType>
    virtual void PseudoSpectral::build( )
    {
      // TODO
    }

  template<ParameterDataType, PhysicsDataType>
    virtual PhysicsDataType& evaluate( 
        ParameterDataType& parameterValues /**< parameter values to evaluate*/
        )
    {
      // TODO
    }

  template<ParameterDataType, PhysicsDataType>
    virtual void refine( )
    {
      m_order = m_order + 1;
    }


  
}
#endif // PSEUDO_SPECTRAL_H


