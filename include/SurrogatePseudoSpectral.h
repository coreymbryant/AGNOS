
#ifndef SURROGATE_PSEUDO_SPECTRAL_H
#define SURROGATE_PSEUDO_SPECTRAL_H
#include "SurrogateModel.h"


namespace AGNOS
{

/********************************************//**
 * \brief SurrogateModel based on Pseudo-spectral projection. 
 *
 * This class provides the framework for constructing surrogate models using
 * non-intrusive spectral projection methods. 
 *
 * As of now only the derived class PseudoSpectralTensorProduct is operational
 * but it could be extended to other methods as well. Both isotropic and
 * non-isotropic polynomial orders are supported. 
 ***********************************************/
  template<class T_S, class T_P>
    class SurrogatePseudoSpectral : public SurrogateModel<T_S,T_P>
  {

    public:

      SurrogatePseudoSpectral( 
          PhysicsFunction<T_S,T_P>&     solutionFunction,
          const std::vector<Parameter*> parameters,
          const unsigned int            order 
          );
      SurrogatePseudoSpectral( 
          PhysicsFunction<T_S,T_P>&         solutionFunction,
          const std::vector<Parameter*>     parameters,
          const std::vector<unsigned int>&  order
          );

      virtual ~SurrogatePseudoSpectral( );

      void build( ) ;
      std::vector<T_P> computeContribution( 
          T_S& integrationPoint, 
          double& integrationWeight,
          std::vector<double>& polyValues
          );
      T_P evaluate( 
          T_S& parameterValues /**< parameter values to evaluate*/
          );

      // Manipulators
      std::vector<unsigned int> getExpansionOrder( ) const;
      unsigned int              getNIntegrationPoints( ) const;
      std::vector<T_S>          getIntegrationPoints( ) const;
      std::vector<double>       getIntegrationWeights( ) const;
      const std::vector<T_P>    getCoefficients( ) const;
      const std::vector< std::vector< unsigned int> > 
                                getIndexSet( ) const;
      std::vector<double>       evaluateBasis( T_S& parameterValues ) const;

      void                      printIntegrationWeights( ) const;
      void                      printIntegrationPoints( ) const;
      void                      printIndexSet( ) const;

    protected:

      std::vector<unsigned int> m_order;  
      unsigned int              m_nIntegrationPoints;
      std::vector<T_S>          m_integrationPoints ;
      std::vector<double>       m_integrationWeights ;

      std::vector<T_P>          m_coefficients;
      std::vector< std::vector<unsigned int> > m_indexSet;



  };


/********************************************//**
 * \brief Constructor for isotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        PhysicsFunction<T_S,T_P>& solutionFunction,
        const std::vector<Parameter*> parameters,
        const unsigned int order
        )
      : SurrogateModel<T_S,T_P>(solutionFunction,parameters)
    {
      m_order = std::vector<unsigned int>(this->m_dimension,order);
    }


/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        PhysicsFunction<T_S,T_P>& solutionFunction,
        const std::vector<Parameter*> parameters,
        const std::vector<unsigned int>& order
        )
      : SurrogateModel<T_S,T_P>(solutionFunction,parameters), 
      m_order(order)
    {
      if (m_order.size() != parameters.size() )
      {
        std::cout 
          << std::endl
          << "\tERROR:"
          << " order vector dimension does not match number of parameters"
          << std::endl
          << std::endl;
        assert(0);
      }
    }

/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::~SurrogatePseudoSpectral( )
    {
    }

/********************************************//**
 * \brief Get the current expansion order
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<unsigned int> SurrogatePseudoSpectral<T_S,T_P>::getExpansionOrder( )
    const
  {
    return m_order;
  }

/********************************************//**
 * \brief Get number of integration points
 ***********************************************/
  template<class T_S, class T_P>
    unsigned int SurrogatePseudoSpectral<T_S,T_P>::getNIntegrationPoints( )
    const
  {
    return m_nIntegrationPoints;
  }

/********************************************//**
 * \brief Get integration weights
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<double> SurrogatePseudoSpectral<T_S,T_P>::getIntegrationWeights( )
    const
  {
    return m_integrationWeights;
  }

/********************************************//**
 * \brief Get the current expansion order
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<T_S> SurrogatePseudoSpectral<T_S,T_P>::getIntegrationPoints( )
    const
  {
    return m_integrationPoints;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P> 
    const std::vector<T_P> 
    SurrogatePseudoSpectral<T_S,T_P>::getCoefficients( ) const
    {
      return m_coefficients;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P> 
    const std::vector< std::vector< unsigned int> > 
    SurrogatePseudoSpectral<T_S,T_P>::getIndexSet( ) const
    {
      return m_indexSet;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P> 
    std::vector<double> SurrogatePseudoSpectral<T_S,T_P>::evaluateBasis( 
        T_S& parameterValue 
        ) const
    {
      unsigned int nTerms = this->m_indexSet.size() ;
      std::vector<double> basisValues( nTerms ,1.);

      for(unsigned int id=0; id < nTerms ; id++)
        for(unsigned int dir=0; dir < this->m_dimension; dir++)
        {
          basisValues[id] 
            *= this->m_parameters[dir]->evalBasisPoly( 
                m_indexSet[id][dir], parameterValue(dir) ) ;
        }

      return basisValues;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogatePseudoSpectral<T_S,T_P>::build( )
    {
      // This is separated from the routine that actually computes contribution
      // so that we can group surrogate models together later and wll that needs
      // to be defined is build routine based on computeContribution( )
      //
      
      std::vector<double> polyValues = evaluateBasis(m_integrationPoints[0]) ;
      // TODO change this to just initialize to one and then add this iteration
      // to loop of integration pts
      m_coefficients = computeContribution( 
            m_integrationPoints[0], 
            m_integrationWeights[0], 
            polyValues
            );
      std::vector<T_P> contrib;

      
      for(unsigned int point=1; point < m_nIntegrationPoints; point++)
      {
        polyValues = evaluateBasis(m_integrationPoints[point]) ;
        contrib = computeContribution( 
            m_integrationPoints[point], 
            m_integrationWeights[point], 
            polyValues
            );
        for(unsigned int coeff=0; coeff < m_coefficients.size(); coeff++)
        {
          m_coefficients[coeff] += contrib[coeff];
        }
      }
      
      return;
    } 

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<T_P> SurrogatePseudoSpectral<T_S,T_P>::computeContribution(
        T_S& integrationPoint, 
        double& integrationWeight,
        std::vector<double>& polyValues
        )
    {
      // TODO pass polyValues or compute locally?
      // TODO this part can be done in parallel 
      T_P solution( this->m_solutionFunction.getImageSize() );
      this->m_solutionFunction.compute( integrationPoint, solution );

      std::vector<T_P> contrib;
      contrib.resize( polyValues.size() ); // polyValues has same size as coeff
      for (unsigned int i=0; i < contrib.size(); i++)
      {
        contrib[i] = solution ;
        for (unsigned int j=0; j < solution.size(); j++)
          contrib[i](j) *= polyValues[i] * integrationWeight ;
      }

      return contrib;
    }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  template<class T_S, class T_P>
    T_P SurrogatePseudoSpectral<T_S,T_P>::evaluate( 
        T_S& parameterValues /**< parameter values to evaluate*/
        )
    {
      T_P surrogateValue = m_coefficients[0];
      std::vector<double> polyValues = evaluateBasis(parameterValues) ;

      // TODO again initialize this somehow and absorb this iteration in loop
      // below
      for(unsigned int comp=0; comp < surrogateValue.size(); comp++)
        surrogateValue(comp) =  m_coefficients[0](comp) * polyValues[0];

      for(unsigned int coeff=1; coeff < m_coefficients.size(); coeff++)
        for(unsigned int comp=0; comp < surrogateValue.size(); comp++)
          surrogateValue(comp) +=  m_coefficients[coeff](comp) * polyValues[coeff];

      return surrogateValue;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
  void SurrogatePseudoSpectral<T_S,T_P>::printIntegrationPoints( ) const
  {
    std::cout << std::endl;
    std::cout << "====================================================" <<
      std::endl;
    std::cout << " Integration points " << std::endl;
    std::cout << "----------------------------------------------------" <<
      std::endl;
    std::cout << "   \\ x        " << std::endl;
    std::cout << "    \\" ;
    for(unsigned int dim=0; dim < this->m_dimension; dim++)
      std::cout << std::setw(12) << "x_" << dim << " " ;
    std::cout << std::endl;
    std::cout << "  id \\  " << std::endl;
    std::cout << "----------------------------------------------------" <<
      std::endl;
    for(int ix=0; ix < m_nIntegrationPoints ; ix++)  
    {
      std::cout << std::setw(5) << ix << " |   ";
      for(int iy=0; iy < this->m_dimension; iy++)
      {
        std::cout << std::scientific << std::setprecision(5) << std::setw(12)
          << m_integrationPoints[ix](iy) << "  ";
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
  void SurrogatePseudoSpectral<T_S,T_P>::printIntegrationWeights( ) const
  {
    double sum = 0.0;
    std::cout << std::endl;
    std::cout << "====================================================" <<
      std::endl;
    std::cout << " Integration weights " << std::endl;
    std::cout << "----------------------------------------------------" <<
      std::endl;
    for(int ip=0; ip < m_nIntegrationPoints; ip++){  
      std::cout << std::setw(5) << ip << "   |   ";
      std::cout << m_integrationWeights[ip] << "  ";
      std::cout << std::endl;
      sum += m_integrationWeights[ip];
    }

    std::cout << std::endl;
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
    void SurrogatePseudoSpectral<T_S,T_P>::printIndexSet( ) const
    {
      std::cout << std::endl;
      std::cout << "====================================================" <<
        std::endl;
      std::cout << " Index Set" << std::endl;
      std::cout << "----------------------------------------------------" <<
      std::endl;
    std::cout << "   \\ dir      " << std::endl;
    std::cout << "    \\        " ;
    for(unsigned int dim=0; dim < this->m_dimension; dim++)
      std::cout << std::setw(4) << "xi_" << dim << " " ;
    std::cout << std::endl;
    std::cout << "  id \\  " << std::endl;
    std::cout << "----------------------------------------------------" <<
      std::endl;
      for (unsigned int i=0; i< m_indexSet.size(); i++)
      {
      std::cout << std::setw(5) << i << "   |   ";
        for (unsigned int j=0; j< m_indexSet[i].size(); j++)
        {
          std::cout << std::setw(5) << m_indexSet[i][j] << " " ;
        }
        std::cout << std::endl;
      }

      return ;
    }
  
}
#endif // SURROGATE_PSEUDO_SPECTRAL_H


