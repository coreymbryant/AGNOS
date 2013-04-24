
#ifndef PSEUDO_SPECTRAL_H
#define PSEUDO_SPECTRAL_H
#include "SurrogateModel.h"
#include "sandia_rules.hpp"

namespace AGNOS
{

  template<class T_S, class T_P>
    class PseudoSpectral : public SurrogateModel<T_S,T_P>
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

      virtual void build( ); 

      virtual T_P& evaluate( 
          T_S& parameterValues /**< parameter values to evaluate*/
          );

      virtual void refine( );

      void setParameterDimension(unsigned int parameterDimension) ;
      unsigned int getParameterDimension( ) const;

    private:

      unsigned int m_order;       // isotropic for now
      unsigned int m_dimension;   // all uniform distributed
      std::vector<double> m_mins; // can have different range though
      std::vector<double> m_maxs;

      std::vector<double> m_coefficients;
      std::vector<double> m_indexSet;

      std::vector< std::vector<double> >  m_quadPoints;
      std::vector<double> m_quadWeights;

  };


  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( )
      : m_order (0), m_dimension(1)
    {
      m_mins.push_back(-1.0);
      m_maxs.push_back(1.0);

      double* points = new double[1];
      double* weights = new double[1];

      webbur::legendre_compute( 1, points, weights );

      m_quadPoints[0].push_back( *points );
      m_quadWeights.push_back(*weights);

      return;
    }

  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( 
        unsigned int order, 
        unsigned int dimension, 
        std::vector<double>& m_mins, 
        std::vector<double>& m_maxs
        )
    {
      PseudoSpectral( );
      return;
    }

  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::~PseudoSpectral()
    {
      return ;
    }


  template<class T_S, class T_P>
    void PseudoSpectral<T_S,T_P>::build( )
    {
      // TODO
    }

  template<class T_S, class T_P>
    T_P& PseudoSpectral<T_S,T_P>::evaluate( 
        T_S& parameterValues /**< parameter values to evaluate*/
        )
    {
      // TODO
    }

  template<class T_S, class T_P>
    void PseudoSpectral<T_S,T_P>::refine( )
    {
      m_order = m_order + 1;
    }


  
}
#endif // PSEUDO_SPECTRAL_H


