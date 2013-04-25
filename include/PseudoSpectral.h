
#ifndef PSEUDO_SPECTRAL_H
#define PSEUDO_SPECTRAL_H
#include "SurrogateModel.h"
#include "sandia_rules.hpp"

namespace AGNOS
{

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    class PseudoSpectral : public SurrogateModel<T_S,T_P>
  {

    public:

      PseudoSpectral( 
          unsigned int order, 
          unsigned int dimension
          );
      PseudoSpectral( 
          unsigned int order, 
          unsigned int dimension, 
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
      unsigned int m_order;       // isotropic for now
      std::vector<double> m_mins; // can have different range though
      std::vector<double> m_maxs;

      std::vector<double> m_coefficients;
      std::vector<double> m_indexSet;

      double**  m_quadPoints;
      double* m_quadWeights;
      unsigned int m_nQuadPoints;

      virtual void constructQuadRule();

      void recurQuad(
          const int dim, const int order, 
          double currentWeights[], double* currentPointsp[] );

  };


/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( )
      : m_order(0), m_dimension(1)
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
        unsigned int order, 
        unsigned int dimension, 
        std::vector<double>& mins, 
        std::vector<double>& maxs
        )
      : m_order(order),m_dimension(dimension), m_mins(mins), m_maxs(maxs)
    {
      constructQuadRule();
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( 
        unsigned int order, 
        unsigned int dimension
        )
      : m_order(order), m_dimension(dimension)
    {
      m_mins.push_back(-1.0);
      m_maxs.push_back(1.0);
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
      m_nQuadPoints = pow(m_order+1,m_dimension);

      m_quadPoints = new double*[m_dimension];
      for(unsigned int i=0; i<m_dimension; i++)
        m_quadPoints[i] = new double[m_nQuadPoints];
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
// TODO adapt to this code
// recursive quad points
/* getQuadPoints() */
/* { */
/*   nQuadPoints_ = pow(order_+1,nDim_) ; */
/*   quadPoints_.Resize( nQuadPoints_, nDim_, NULL); */

/*   UQTKArray1D<double> points1D = p_basis_->quadPoints_; */
  
/*   for(int dit=0; dit < nDim_; dit++) */
/*   { */
/*     int total = 0; */
/*     while(total < nQuadPoints_) */
/*     { */
/*       for(int pit=0; pit < points1D.XSize(); pit++) */
/*         for(int rit=0; rit < pow(order_+1,nDim_-(dit+1)); rit++) */
/*         { */
/*           quadPoints_(total,dit)   = points1D(pit); */
/*           if (UQTKVERBOSE) */
/*             cout << "qPoint ( " << total << " , " << dit << " ) = " << */
/*               quadPoints_(total,dit) << endl; */
/*           total++; */
/*         } */
/*       } */
/*   } */
/*   return quadPoints_; */
/* } */
  
/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void PseudoSpectral<T_S,T_P>::recurQuad(
        const int dim, const int order, 
        double currentWeights[], double* currentPoints[] )
    {
      double* oneDimQuadWeights = new double[order+1];
      double* oneDimQuadPoints = new double[order+1];
      webbur::legendre_compute( order+1, oneDimQuadPoints, oneDimQuadWeights );

      double scaling = (m_maxs[dim-1]-m_mins[dim-1])/2.0;
      double midpoint = (m_maxs[dim-1]+m_mins[dim-1])/2.0;
          std::cout << scaling << std::endl;

      if (dim == 1) 
      {

        for(unsigned int id=0; id < order+1; id++)
        {
          /* currentPoints[dim-1][id] = midpoint */ 
          /*   + oneDimQuadPoints[id] * scaling; */
            
          currentWeights[id] = oneDimQuadWeights[id] * scaling;
          std::cout << currentWeights[id] << std::endl;
          std::cout << scaling << std::endl;
        }
      }

      else 
      {
        recurQuad(dim-1,order,currentWeights,currentPoints);
        
        // keep track of how many array elements are non-zero
        int prevSize = pow(order+1,dim-1) ;

        for(int outer=0; outer < prevSize; outer++)
        {
          currentWeights[(outer+1)*(order+1)] = currentWeights[]
          for(int inner=0; inner < order+1 ; inner++)
          {
            currentWeights[(outer)*(order+1) + inner] = 
              currentWeights[outer] * oneDimQuadWeights[inner] * scaling;
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
            << m_quadPoints[iy][ix] << "  ";
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
      }
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
      m_order = m_order + 1;
    }


  
}
#endif // PSEUDO_SPECTRAL_H


