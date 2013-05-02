// TODO need to change double** into some sort of array of T_P

#ifndef PSEUDO_SPECTRAL_H
#define PSEUDO_SPECTRAL_H
#include "SurrogateModel.h"


namespace AGNOS
{

/********************************************//**
 * \brief Pseudo-spectral projection surrogate model
 *
 * This class provides the framework for constructing surrogate models using
 * non-intrusive spectral projection methods. 
 *
 * As of now it is only Tensor product quadrature to comptue coefficients but
 * could be extended to other methods as well. Both isotropic and non-isotropic
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
          PhysicsFunction<T_S,T_P>&     solutionFunction,
          const std::vector<Parameter*> parameters,
          const unsigned int            order 
          );
      PseudoSpectral( 
          PhysicsFunction<T_S,T_P>&         solutionFunction,
          const std::vector<Parameter*>     parameters,
          const std::vector<unsigned int>&  order
          );

      virtual ~PseudoSpectral( );

      void build( ) ;
      T_P pointContribution( T_S& integrationPoint );
      T_P evaluate( 
          T_S& parameterValues /**< parameter values to evaluate*/
          );

      // Manipulators
      std::vector<unsigned int> getExpansionOrder( ) const;
      unsigned int              getNIntegrationPoints( ) const;
      std::vector<T_S>          getIntegrationPoints( ) const;
      std::vector<double>       getIntegrationWeights( ) const;
      std::vector< std::vector<double> > getBasisEvaluations( ) const;

      void                      printIntegrationWeights( ) const;
      void                      printIntegrationPoints( ) const;

    protected:

      std::vector<unsigned int> m_order;  
      unsigned int              m_nIntegrationPoints;
      std::vector<T_S>          m_integrationPoints ;
      std::vector<double>       m_integrationWeights ;
      std::vector< std::vector<double> >       m_basisEvaluations ;

      // TODO 
      std::vector<T_P>          m_coefficients;
      std::vector< std::vector<unsigned int> > m_indexSet;



  };


/********************************************//**
 * \brief Constructor for isotropic order
 ***********************************************/
  template<class T_S, class T_P>
    PseudoSpectral<T_S,T_P>::PseudoSpectral( 
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
    PseudoSpectral<T_S,T_P>::PseudoSpectral( 
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
    PseudoSpectral<T_S,T_P>::~PseudoSpectral( )
    {
    }

/********************************************//**
 * \brief Get the current expansion order
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<unsigned int> PseudoSpectral<T_S,T_P>::getExpansionOrder( )
    const
  {
    return m_order;
  }

/********************************************//**
 * \brief Get number of integration points
 ***********************************************/
  template<class T_S, class T_P>
    unsigned int PseudoSpectral<T_S,T_P>::getNIntegrationPoints( )
    const
  {
    return m_nIntegrationPoints;
  }

/********************************************//**
 * \brief Get integration weights
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<double> PseudoSpectral<T_S,T_P>::getIntegrationWeights( )
    const
  {
    return m_integrationWeights;
  }

/********************************************//**
 * \brief Get the current expansion order
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<T_S> PseudoSpectral<T_S,T_P>::getIntegrationPoints( )
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
    void PseudoSpectral<T_S,T_P>::build( )
    {
      // This is separated from the routine that actually computes contribution
      // so that we can group surrogate models together later and wll that needs
      // to be defined is build routine based on pointContribution( )

      // m_integrationWeights
      // poly value

      for(unsigned int point=0; point < m_nIntegrationPoints; point++)
        pointContribution( 
            m_integrationWeights[point], 
            m_integrationPoints[point], 
            polyValue 
            );
    } 

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    void PseudoSpectral<T_S,T_P>::pointContribution(
        T_S& integrationPoint
        )
    {

      for (unsigned int j=0; j < m_coefficients.size(); j++)
      {
        T_P tmpSolution;
        m_solutionFunction.compute( integrationPoint, tmpSolution );

        std::vector<double> 
        if (m_coefficients[j].size() == 0)
          m_coefficients[j] = tmpSolution  ;
        else
          m_coefficients[j] += 

      }
      /* T_P solVec; */ 

      /* unsigned int nQuadPoints = this->m_quadRule->getNQuadPoints(); */
      /* std::vector<T_S> parameterValues */
      /*   = this->m_quadRule->getQuadPoints(); */

      /* for (unsigned int j=0; j < nQuadPoints; j++) */
      /* { */
      /*   this->m_solutionFunction.compute(paramValues[j],solVec); */

      /*   for(unsigned int coeff=0; coeff < nQuadPoints; coeff++) */
      /*   { */
      /*     std::cout << "solVec[0] = " << solVec[0] << std::endl; */
      /*     std::cout << "poly = " */ 
      /*       << (this->m_parameters[0]->evalBasisPoly(coeff,paramVec[0]) ) */ 
      /*       << std::endl; */
      /*     std::cout << "weight = " */ 
      /*       << (this->m_quadRule->getQuadWeights())[0] */
      /*       << std::endl; */

      /*     for(unsigned int solId; solId < numberOfSolutionVectors; solId++) */
      /*       m_coefficients[coeff][solId][0] +=  solVec[solId][0] */
      /*         * (this->m_parameters[0]->evalBasisPoly(coeff,paramVec[0][0]) ) */ 
      /*         * (this->m_quadRule->getQuadWeights())[0] ; */
      /*   } */
      /*   std::cout << "coeff[" << j << "][0] = " << m_coefficients[j][0][0] << std::endl; */
      /* } */
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S, class T_P>
    T_P PseudoSpectral<T_S,T_P>::evaluate( 
        T_S& parameterValues /**< parameter values to evaluate*/
        )
    {
      /* T_P currSol; */
      /* // TODO */
      /* for(unsigned int coeff=0; coeff < m_nQuadPoints; coeff++) */
      /*   /1* currSol += coeff * Poly ; *1/ */
      /*   currSol = 0; */
      T_P dummy;

      return dummy;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
  void PseudoSpectral<T_S,T_P>::printIntegrationPoints( ) const
  {
    std::cout << std::endl;
    std::cout << "====================================================" <<
      std::endl;
    std::cout << " Integration points " << std::endl;
    std::cout << "----------------------------------------------------" <<
      std::endl;
    std::cout << "  \\ x        " << std::endl;
    std::cout << "   \\" ;
    for(unsigned int dim=0; dim < this->m_dimension; dim++)
      std::cout << std::setw(12) << "x_" << dim << " " ;
    std::cout << std::endl;
    std::cout << " id \\  " << std::endl;
    std::cout << "----------------------------------------------------" <<
      std::endl;
    for(int ix=0; ix < m_nIntegrationPoints ; ix++)  
    {
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
  void PseudoSpectral<T_S,T_P>::printIntegrationWeights( ) const
  {
    double sum = 0.0;
    std::cout << std::endl;
    std::cout << "====================================================" <<
      std::endl;
    std::cout << " Integration weights " << std::endl;
    std::cout << "----------------------------------------------------" <<
      std::endl;
    for(int ip=0; ip < m_nIntegrationPoints; ip++){  
      std::cout << ip << "   | ";
      std::cout << m_integrationWeights[ip] << "  ";
      std::cout << std::endl;
      sum += m_integrationWeights[ip];
    }
    std::cout << "Sum = " << sum << std::endl;
    std::cout << std::endl;
    return;
  }
  
}
#endif // PSEUDO_SPECTRAL_H


