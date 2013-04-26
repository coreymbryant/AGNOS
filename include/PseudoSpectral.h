
#ifndef PSEUDO_SPECTRAL_H
#define PSEUDO_SPECTRAL_H
#include "SurrogateModel.h"
#include "sandia_rules.hpp"

namespace AGNOS
{

/********************************************//**
 * \brief Pseudo-spectral projection surrogate model
 *
 * This class provides the framework for constructing surrogate models using
 * non-intrusive spectral projection methods. Both isotropic and non-isotropic
 * polynomial orders are supported. 
 *
 * Only Uniform Random variables at the moment
 * 
 ***********************************************/
  template<class T_S, class T_P>
    class PseudoSpectral : public SurrogateModel<T_S,T_P>
  {

    public:

      PseudoSpectral( 
          unsigned int dimension,
          unsigned int order 
          );
      PseudoSpectral( 
          unsigned int dimension,
          std::vector<unsigned int> order
          );
      PseudoSpectral( 
          unsigned int dimension, 
          unsigned int order, 
          std::vector<double>& m_mins, 
          std::vector<double>& m_maxs
          );
      PseudoSpectral( 
          unsigned int dimension, 
          std::vector<unsigned int> order, 
          std::vector<double>& m_mins, 
          std::vector<double>& m_maxs
          );

      PseudoSpectral( );

      ~PseudoSpectral( );

      void build( ); 

      T_P& evaluate( 
          T_S& parameterValues /**< parameter values to evaluate*/
          );

      void refine( );

      void setParameterDimension(unsigned int parameterDimension) ;
      unsigned int getParameterDimension( ) const;

      void printQuadWeights( ) const;
      void printQuadPoints( ) const;

    private:

      unsigned int m_dimension;   // all uniform distributed
      std::vector<unsigned int> m_order; // anisotropic
      std::vector<double> m_mins; // can have different range though
      std::vector<double> m_maxs;

      std::vector<double> m_coefficients;
      std::vector<double> m_indexSet;

      double**  m_quadPoints;
      double* m_quadWeights;
      unsigned int m_nQuadPoints;

      void constructQuadRule();

      // We can use this to derive other classes for non-uniform RV (if needed)
      virtual void oneDimQuadRule(
        unsigned int order, 
        double oneDimQuadPoints[], double oneDimQuadWeights[] );

      void recurQuad(
          const int dim, const std::vector<unsigned int> order, 
          double currentWeights[], double* currentPoints[] );

  };


/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( )
      : m_dimension(1)
    {
      m_order.push_back(0);
      m_mins.push_back(-1.0);
      m_maxs.push_back(1.0);
      constructQuadRule();
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( 
        unsigned int dimension,
        unsigned int order
        )
      : m_dimension(dimension)
    {
      m_order = std::vector<unsigned int>(dimension,order);
      m_mins.push_back(-1.0);
      m_maxs.push_back(1.0);
      constructQuadRule();
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( 
        unsigned int dimension,
        std::vector<unsigned int> order
        )
      : m_order(order), m_dimension(dimension)
    {
      m_mins.push_back(-1.0);
      m_maxs.push_back(1.0);
      constructQuadRule();
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( 
        unsigned int dimension, 
        unsigned int order, 
        std::vector<double>& mins, 
        std::vector<double>& maxs
        )
      : m_dimension(dimension), m_mins(mins), m_maxs(maxs)
    {
      m_order = std::vector<unsigned int>(dimension,order);
      constructQuadRule();
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( 
        unsigned int dimension, 
        std::vector<unsigned int> order, 
        std::vector<double>& mins, 
        std::vector<double>& maxs
        )
      : m_order(order),m_dimension(dimension), m_mins(mins), m_maxs(maxs)
    {
      constructQuadRule();
    }

  
/********************************************//**
 * \brief Rountine to populate quadPoints and quadWeights
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
    void PseudoSpectral<T_S,T_P>::constructQuadRule( )
    {
      m_nQuadPoints = m_order[0]+1;
      for(unsigned int i=1; i< m_dimension; i++)
        m_nQuadPoints *= (m_order[i]+1);

      m_quadPoints = new double*[m_nQuadPoints];
      for(unsigned int i=0; i<m_nQuadPoints; i++)
        m_quadPoints[i] = new double[m_dimension];
      m_quadWeights = new double[m_nQuadPoints];

      recurQuad(m_dimension,m_order,m_quadWeights,m_quadPoints);
      /* setQuadPoints( ); */


      return ;
    }
  
/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void PseudoSpectral<T_S,T_P>::recurQuad(
        const int dim, const std::vector<unsigned int> order, 
        double currentWeights[], double* currentPoints[] )
    {
      double* oneDimQuadWeights = new double[order[dim-1]+1];
      double* oneDimQuadPoints = new double[order[dim-1]+1];
      oneDimQuadRule( order[dim-1]+1, oneDimQuadPoints, oneDimQuadWeights );

      double scaling = (m_maxs[dim-1]-m_mins[dim-1])/2.0;
      double midpoint = (m_maxs[dim-1]+m_mins[dim-1])/2.0;

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
        recurQuad(dim-1,order,currentWeights,currentPoints);
        
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
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void PseudoSpectral<T_S,T_P>::oneDimQuadRule(
        unsigned int order, 
        double oneDimQuadPoints[], double oneDimQuadWeights[] )
    {
      webbur::legendre_compute( order, oneDimQuadPoints, oneDimQuadWeights );
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::~PseudoSpectral()
    {
      delete [] m_quadPoints;
      delete [] m_quadWeights;
    }


/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void PseudoSpectral<T_S,T_P>::printQuadPoints( ) const
    {
      std::cout << std::endl;
      std::cout << "====================================================" <<
        std::endl;
      std::cout << " Quadrature points " << std::endl;
      std::cout << "----------------------------------------------------" <<
        std::endl;
      std::cout << "  \\ x        " << std::endl;
      std::cout << "   \\" ;
      for(unsigned int dim=0; dim < m_dimension; dim++)
        std::cout << std::setw(12) << "x_" << dim << " " ;
      std::cout << std::endl;
      std::cout << " id \\  " << std::endl;
      std::cout << "----------------------------------------------------" <<
        std::endl;
      for(int ix=0; ix < m_nQuadPoints ; ix++)  
      {
        std::cout << std::setw(3) << ix << "  |  " ;
        for(int iy=0; iy < m_dimension; iy++)
        {
          std::cout << std::scientific << std::setprecision(5) << std::setw(12)
            << m_quadPoints[ix][iy] << "  ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
      return;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
    void PseudoSpectral<T_S,T_P>::printQuadWeights( ) const
    {
      double sum = 0.0;
      std::cout << std::endl;
      std::cout << "====================================================" <<
        std::endl;
      std::cout << " Quadrature weights " << std::endl;
      std::cout << "----------------------------------------------------" <<
        std::endl;
      for(int ip=0; ip < m_nQuadPoints; ip++){  
        std::cout << ip << "   | ";
        std::cout << m_quadWeights[ip] << "  ";
        std::cout << std::endl;
        sum += m_quadWeights[ip];
      }
      std::cout << "Sum = " << sum << std::endl;
      std::cout << std::endl;
      return;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void PseudoSpectral<T_S,T_P>::build( )
    {
      // TODO
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    T_P& PseudoSpectral<T_S,T_P>::evaluate( 
        T_S& parameterValues /**< parameter values to evaluate*/
        )
    {
      // TODO
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void PseudoSpectral<T_S,T_P>::refine( )
    {
      for(unsigned int i=0; i<m_dimension; i++)
        m_order[i]++;
    }


  
}
#endif // PSEUDO_SPECTRAL_H


