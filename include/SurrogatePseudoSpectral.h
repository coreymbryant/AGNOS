
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

      // single physics function constructors
      SurrogatePseudoSpectral( 
          PhysicsFunction<T_S,T_P>*     solutionFunction,
          const std::vector<Parameter*> parameters,
          const unsigned int            order 
          );
      SurrogatePseudoSpectral( 
          PhysicsFunction<T_S,T_P>*         solutionFunction,
          const std::vector<Parameter*>     parameters,
          const std::vector<unsigned int>&  order
          );

      // multiple physics function constructors
      SurrogatePseudoSpectral( 
          std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
          const std::vector<Parameter*>                       parameters,
          const unsigned int                                  order 
          );
      SurrogatePseudoSpectral( 
          std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
          const std::vector<Parameter*>                       parameters,
          const std::vector<unsigned int>&                    order
          );

      virtual ~SurrogatePseudoSpectral( );

      void build( ) ;
      std::map< std::string, std::vector<T_P> > computeContribution( 
          std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
          T_S&                                                integrationPoint, 
          double&                                             integrationWeight,
          std::vector<double>                                 polyValues
          );

      using SurrogateModel<T_S,T_P>::evaluate; 
      std::map<std::string, T_P> evaluate( 
          std::vector< std::string >  solutionNames,
          T_S&                        parameterValues 
          );

      // Manipulators
      unsigned int              getNIntegrationPoints( ) const;
      std::vector<T_S>          getIntegrationPoints( ) const;
      std::vector<double>       getIntegrationWeights( ) const;
      const std::vector< std::vector< unsigned int> > 
                                getIndexSet( ) const;
      std::vector<double>       evaluateBasis( T_S& parameterValues ) const;
      // TODO change these to 'pretty print' or something
      void                      printIntegrationWeights( ) const;
      void                      printIntegrationPoints( ) const;
      void                      printIndexSet( ) const;

    protected:

      unsigned int              m_nIntegrationPoints;
      std::vector<T_S>          m_integrationPoints ;
      std::vector<double>       m_integrationWeights ;

      std::vector< std::vector<unsigned int> > m_indexSet;



  };


/********************************************//**
 * \brief Constructor for isotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        PhysicsFunction<T_S,T_P>* solutionFunction,
        const std::vector<Parameter*> parameters,
        const unsigned int order
        )
      : SurrogateModel<T_S,T_P>(solutionFunction,parameters,order)
    {
    }


/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        PhysicsFunction<T_S,T_P>* solutionFunction,
        const std::vector<Parameter*> parameters,
        const std::vector<unsigned int>& order
        )
      : SurrogateModel<T_S,T_P>(solutionFunction,parameters,order) 
    {
    }

/********************************************//**
 * \brief Constructor for isotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
        const std::vector<Parameter*>                       parameters,
        const unsigned int                                  order
        )
      : SurrogateModel<T_S,T_P>(solutionFunction,parameters,order)
    {
    }


/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
        const std::vector<Parameter*>                       parameters,
        const std::vector<unsigned int>&                    order
        )
      : SurrogateModel<T_S,T_P>(solutionFunction,parameters,order) 
    {
    }

/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::~SurrogatePseudoSpectral( )
    {
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
      
      
      // INITIAL INTEGRATION POINT TO SET SIZES 
      
      // TODO would be nice if we could initialize empty vector of correct size
      // and then run all integration pts in parallel. This would require some
      // manipulator to return physics vector size. Instead we can just make one
      // call to computeContribution to set sizes and then run the rest of the
      // jobs in parrael. I guess we could probably do all computations in
      // parallel and just deal with the the size issue before combining back to
      // one. 
      std::vector<double> polyValues;
      typename std::map< std::string, std::vector<T_P> >::iterator id;
      
      for(unsigned int point=0; point < m_nIntegrationPoints; point++)
      {
        polyValues = evaluateBasis(m_integrationPoints[point]) ;
        std::map< std::string, std::vector<T_P> > contrib = computeContribution( 
            this->m_solutionFunction,
            m_integrationPoints[point], 
            m_integrationWeights[point], 
            polyValues
            );

        // if this is the first integration point initialize coeff to its value
        if (point==0)
          this->m_coefficients = contrib;          

        else
          // TODO change this
          for(id=this->m_coefficients.begin(); id!=this->m_coefficients.end(); id++)
            for(unsigned int coeff=0; coeff < (id->second).size(); coeff++)
              (id->second)[coeff] += (contrib[id->first])[coeff];
      }
      
      return;
    } 

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    std::map< std::string, std::vector<T_P> >
    SurrogatePseudoSpectral<T_S,T_P>::computeContribution(
        std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
        T_S&                                                integrationPoint, 
        double&                                             integrationWeight,
        std::vector<double>                                 polyValues
        )
    {
      // TODO this part can be done in parallel 
      
      // I think its better to pass vector of polyValues than a pointer to
      // SurrogateModel object, since in parallel their would be multiple nodes
      // operating on same object (even though they would just be referencing
      // not altering that object?? - yes pointer would be to object stored in
      // master node not working node)

      std::map< std::string, std::vector<T_P> > contrib;
      T_P solution ;

      
      typename std::map< std::string, PhysicsFunction<T_S,T_P>* >::iterator id;
      for (id=solutionFunction.begin(); id!=solutionFunction.end(); id++)
      {
        
        // get solution for this integration point
        id->second->compute( integrationPoint, solution );

        // compute contribution of current solution to overall coeff vector
        std::vector<T_P> dummyVec(polyValues.size(), solution);
        contrib.insert( std::pair< std::string, std::vector<T_P> >(
              id->first, dummyVec  ) ); 
              /* id->first, std::vector<T_P>(polyValues.size(),solution)  ) ); */ 

        for (unsigned int i=0; i < contrib[id->first].size(); i++)
          for (unsigned int j=0; j < contrib[id->first][i].size(); j++)
          {
            contrib[id->first][i](j) *= polyValues[i] * integrationWeight ;
          }

        /* contrib[nSol].resize( polyValues.size() ); // polyValues has same size as coeff */
        /* for (unsigned int i=0; i < contrib[nSol].size(); i++) */
        /* { */
        /*   contrib[nSol][i] = solution[nSol] ; */
        /*   for (unsigned int j=0; j < solution[nSol].size(); j++) */
        /*     contrib[nSol][i](j) *= polyValues[i] * integrationWeight ; */
        /* } */
      }



      return contrib;
    }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  template<class T_S, class T_P>
    std::map<std::string, T_P> SurrogatePseudoSpectral<T_S,T_P>::evaluate( 
        std::vector< std::string > solutionNames,
        T_S& parameterValues /**< parameter values to evaluate*/
        )
    {
      std::vector<double> polyValues = evaluateBasis(parameterValues) ;

      // TODO again initialize this somehow and absorb this iteration in loop
      // below
      std::map< std::string, T_P> surrogateValue;

      for (unsigned int i=0; i < solutionNames.size(); i++)
      {
        std::string id = solutionNames[i];
        surrogateValue.insert( 
            std::pair< std::string, T_P>( id, (this->m_coefficients[id])[0] ) 
            );
        for(unsigned int comp=0; comp < (surrogateValue[id]).size(); comp++)
          (surrogateValue[id])(comp) 
            = (this->m_coefficients[id])[0](comp) * polyValues[0];

        for(unsigned int coeff=1; coeff < (this->m_coefficients[id]).size(); coeff++)
          for(unsigned int comp=0; comp < (surrogateValue[id]).size(); comp++)
            (surrogateValue[id])(comp) 
              +=  (this->m_coefficients[id])[coeff](comp) * polyValues[coeff];
      }

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


