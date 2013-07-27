
#ifndef SURROGATE_PSEUDO_SPECTRAL_H
#define SURROGATE_PSEUDO_SPECTRAL_H
#include "SurrogateModel.h"
#include "QuadratureTensorProduct.h"


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
          const Communicator*               comm,
          PhysicsFunction<T_S,T_P>*     solutionFunction,
          const std::vector<Parameter*> parameters,
          const unsigned int            order 
          );
      SurrogatePseudoSpectral( 
          const Communicator*               comm,
          PhysicsFunction<T_S,T_P>*         solutionFunction,
          const std::vector<Parameter*>     parameters,
          const std::vector<unsigned int>&  order
          );

      // multiple physics function constructors
      SurrogatePseudoSpectral( 
          const Communicator*               comm,
          std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
          const std::vector<Parameter*>                       parameters,
          const unsigned int                                  order 
          );
      SurrogatePseudoSpectral( 
          const Communicator*               comm,
          std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
          const std::vector<Parameter*>                       parameters,
          const std::vector<unsigned int>&                    order
          );

      virtual ~SurrogatePseudoSpectral( );

      void build( ) ;
      std::map< std::string, T_P > computeContribution( 
          std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
          T_S&                                                integrationPoint, 
          double&                                             integrationWeight
          );

      using SurrogateModel<T_S,T_P>::evaluate; 
      std::map<std::string, T_P> evaluate( 
          std::vector< std::string >  solutionNames,
          T_S&                        parameterValues 
          );

      using SurrogateModel<T_S,T_P>::l2Norm;
      std::map< std::string, T_P> l2Norm(
        std::vector< std::string > solutionNames
        );
      double l2NormDifference(
          SurrogateModel<T_S,T_P>& comparisonModel,
          std::string solutionName );

      // Manipulators
      unsigned int              getNIntegrationPoints( ) const;
      std::vector<T_S>          getIntegrationPoints( ) const;
      std::vector<double>       getIntegrationWeights( ) const;
      const std::vector< std::vector< unsigned int> > 
                                getIndexSet( ) const;
      std::vector<double>       evaluateBasis( 
          std::vector< std::vector<unsigned int> > indexSet,
          T_S& parameterValues ) const;

      void                      printIntegrationWeights( std::ostream& out ) const;
      void                      printIntegrationPoints( std::ostream& out ) const;
      void                      printIndexSet( std::ostream& out ) const;
      void                      prettyPrintIntegrationWeights( ) const;
      void                      prettyPrintIntegrationPoints( ) const;
      void                      prettyPrintIndexSet( ) const;



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
        const Communicator*               comm,
        PhysicsFunction<T_S,T_P>* solutionFunction,
        const std::vector<Parameter*> parameters,
        const unsigned int order
        )
      : SurrogateModel<T_S,T_P>(comm,solutionFunction,parameters,order) 
    {
    }


/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        const Communicator*               comm,
        PhysicsFunction<T_S,T_P>* solutionFunction,
        const std::vector<Parameter*> parameters,
        const std::vector<unsigned int>& order
        )
      : SurrogateModel<T_S,T_P>(comm,solutionFunction,parameters,order) 
    {
    }

/********************************************//**
 * \brief Constructor for isotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        const Communicator*               comm,
        std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
        const std::vector<Parameter*>                       parameters,
        const unsigned int                                  order
        )
      : SurrogateModel<T_S,T_P>(comm,solutionFunction,parameters,order) 
    {
    }


/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        const Communicator*               comm,
        std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
        const std::vector<Parameter*>                       parameters,
        const std::vector<unsigned int>&                    order
        )
      : SurrogateModel<T_S,T_P>(comm,solutionFunction,parameters,order) 
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
        std::vector< std::vector< unsigned int > > indexSet,
        T_S& parameterValue 
        ) const
    {
      unsigned int nTerms = indexSet.size() ;
      std::vector<double> basisValues( nTerms ,1.);

      for(unsigned int id=0; id < nTerms ; id++)
        for(unsigned int dir=0; dir < this->m_dimension; dir++)
        {
          basisValues[id] 
            *= this->m_parameters[dir]->evalBasisPoly( 
                indexSet[id][dir], parameterValue(dir) ) ;
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
      unsigned int totalNCoeff = this->m_totalNCoeff;
      std::vector< std::vector<double> > polyValues;
      polyValues.reserve(totalNCoeff);
      typename std::map< std::string, PhysicsFunction<T_S,T_P>* >::iterator id;
      
      // output progress
      if( this->m_comm->rank() == 0)
      {
        std::cout << std::endl;
        std::cout << " Starting build of surrogate model for: { " ;
        for (id=this->m_solutionFunction.begin();
            id!=this->m_solutionFunction.end(); id++)
        {
          if (id!=this->m_solutionFunction.begin())
            std::cout << " , " ;
          std::cout << id->first;
        }
        std::cout << " }" << std::endl;
      }
      
      // Int points and coeff are currently separated in case we want to
      // implememnt higher order integration rule in the future
      // assign the appropriate int points to each processor
      unsigned int nPts  
        =  m_nIntegrationPoints / this->m_comm->size()  
        + ( this->m_comm->rank() < (m_nIntegrationPoints % this->m_comm->size() ) ) ;
      unsigned int intPtsStart = std::min(this->m_comm->rank(), m_nIntegrationPoints-1);



      // assign the appropriate coeff to each processor
      unsigned int nCoeffs  
        =  totalNCoeff / this->m_comm->size()  
        + ( this->m_comm->rank() < (totalNCoeff % this->m_comm->size() ) ) ;
      unsigned int coeffStart = std::min(this->m_comm->rank(), totalNCoeff-1);


      // initiaize contribution vector
      std::vector< std::map< std::string, T_P > > myContribs;
      myContribs.reserve(nPts);

      this->m_comm->barrier();
      // compute necessary info for each solution function
      for (id=this->m_solutionFunction.begin();
          id!=this->m_solutionFunction.end(); id++)
      {
        id->second->computeData( m_integrationPoints );
      }
      this->m_comm->barrier();


      std::cout << "order:" << this->m_order[0] << std::endl;
      if( this->m_comm->rank() == 0)
        std::cout << "     --> Solving at " << m_nIntegrationPoints 
          << " integration points " << std::endl;

      // solve for my integration points
      for(unsigned int pt=0; pt < nPts; pt++)
      {
        /* std::cout << "test: beginning of pt" << std::endl; */

        // clear out old solutions
        for (id=this->m_solutionFunction.begin();
            id!=this->m_solutionFunction.end(); id++)
        {
          // set data for each solution function at this integration pt
          id->second->setData( intPtsStart + pt*this->m_comm->size() );
        }


        /* std::cout << "test: pt (pre computeContribution)" << std::endl; */
        // compute the contribution
        myContribs.push_back( computeContribution( 
            this->m_solutionFunction,
            m_integrationPoints[intPtsStart + pt*this->m_comm->size() ], 
            m_integrationWeights[intPtsStart + pt*this->m_comm->size() ]
            ) );
        /* std::cout << "test: pt (post computeContribution)" << std::endl; */
        polyValues.push_back(
            evaluateBasis( this->m_indexSet, m_integrationPoints[intPtsStart +
              pt*this->m_comm->size()]) 
            );
        /* std::cout << "test: end of pt" << std::endl; */
      }

      // need to know size of solution vector on all processes (in case some
      // don't have any work but will still be called in reduce operation)
      unsigned int tempSize;
      this->m_solSize.clear();
      for (id=this->m_solutionFunction.begin();
          id!=this->m_solutionFunction.end(); id++)
      {
        if (this->m_comm->rank() == 0)
          tempSize = myContribs[0][id->first].size();

        this->m_comm->broadcast(tempSize);

        this->m_solSize.insert( std::pair<std::string,unsigned int>(
              id->first,tempSize ) );

      }



      // WAIT FOR ALL PROCESSES TO CATCH UP
      this->m_comm->barrier();
      if( this->m_comm->rank() == 0)
        std::cout << "     --> Computing " << totalNCoeff << " coefficients" <<
          std::endl;


      // --------------
      // compute my contribution to each coefficient
      // sum all processes contributions
      // save if its one of my coefficients
      // --------------
      
      // clear old coefficients
      this->m_coefficients.clear();

      //-- loop through sols
      for (id=this->m_solutionFunction.begin();
          id!=this->m_solutionFunction.end(); id++)
      {
        // temporary storage for my coefficients for this sol
        std::vector<T_P> solCoefficientVectors;
        solCoefficientVectors.reserve(nCoeffs);

        for (unsigned int c=0; c < totalNCoeff; c++)
        {
          std::vector<double> sumContrib(this->m_solSize[id->first],0.);
          
          // compute current rank's contribution
          for(unsigned int pt=0; pt < nPts; pt++)
            for(unsigned int i=0; i <  sumContrib.size(); i++)
              sumContrib[i] += myContribs[pt][id->first](i) *
                polyValues[pt][c];

          // sum all contributions
          this->m_comm->sum( sumContrib );

          // store if its one of mine
          if ( (c % this->m_comm->size()) == this->m_comm->rank() )
            solCoefficientVectors.push_back( T_P(sumContrib) ) ;

        } // c

        if (solCoefficientVectors.size() > 0)
        {
          this->m_coefficients.insert(
              std::pair< std::string, std::vector<T_P> >(id->first,
                solCoefficientVectors) ) ;
        }

      } // id

      return;
    } 

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    std::map< std::string, T_P >
    SurrogatePseudoSpectral<T_S,T_P>::computeContribution(
        std::map< std::string, PhysicsFunction<T_S,T_P>* >  solutionFunction,
        T_S&                                                integrationPoint, 
        double&                                             integrationWeight
        )
    {

      std::cout << "test: computeContribution( ) beginning" << std::endl;
      

      std::map< std::string, T_P > contrib;
      T_P solution ;

      
      typename std::map< std::string, PhysicsFunction<T_S,T_P>* >::iterator id;
      for (id=solutionFunction.begin(); id!=solutionFunction.end(); id++)
      {
        
        std::cout << "test: computeContribution( ) pre compute()" 
          << id->first << std::endl;
        // get solution for this integration point
        id->second->compute( integrationPoint, solution );
        std::cout << "test: computeContribution( ) post compute()" << std::endl;
        solution.scale( integrationWeight );
        contrib.insert( std::pair< std::string, T_P >(
              id->first, solution  ) ); 

      }




      std::cout << "test: computeContribution( ) ending" << std::endl;
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
      unsigned int totalNCoeff = this->m_totalNCoeff;
      std::map< std::string, T_P> surrogateValue;

      // get polyValues
      std::vector<double> polyValues =
        evaluateBasis(this->m_indexSet,parameterValues) ;

      // get starting coeff for current rank
      unsigned int coeffStart = std::min(this->m_comm->rank(), totalNCoeff-1);

      // loop over all solutions requested
      for (unsigned int n=0; n < solutionNames.size(); n++)
      {
        std::string id = solutionNames[n];

        std::vector<double> sumContrib(this->m_solSize[id],0.);

        for(unsigned int c=0; c < this->m_coefficients[id].size(); c++)
        {
          // compute current rank's contribution
          for(unsigned int i=0; i < this->m_coefficients[id][c].size(); i++)
            sumContrib[i] += this->m_coefficients[id][c](i) *
              polyValues[coeffStart + c*(this->m_comm->size())];

        } // coefficients

        // sum all contributions
        this->m_comm->sum( sumContrib );

        surrogateValue.insert(
            std::pair< std::string, T_P >(id, T_P(sumContrib) ) 
            ) ;

      } // solNames

      return surrogateValue;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
  void SurrogatePseudoSpectral<T_S,T_P>::printIntegrationPoints( 
      std::ostream& out ) const
  {
    if (this->m_comm->rank() == 0)
    {
      out << std::endl;
      out << "#====================================================" <<
        std::endl;
      out << "# Integration points " << std::endl;
      out << "#----------------------------------------------------" <<
        std::endl;
      for(int ix=0; ix < m_nIntegrationPoints ; ix++)  
      {
        for(int iy=0; iy < this->m_dimension; iy++)
        {
          out << std::scientific << std::setprecision(5) << std::setw(12)
            << m_integrationPoints[ix](iy) << "  ";
        }
        out << std::endl;
      }
      out << "#----------------------------------------------------" <<
        std::endl;
    }

    return;
  }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
  void SurrogatePseudoSpectral<T_S,T_P>::prettyPrintIntegrationPoints( ) const
  {
    if (this->m_comm->rank() == 0)
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
    }
    return;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
  void SurrogatePseudoSpectral<T_S,T_P>::printIntegrationWeights( 
      std::ostream& out ) const
  {
    if (this->m_comm->rank() == 0)
    {
      double sum = 0.0;
      out << std::endl;
      out << "#====================================================" <<
        std::endl;
      out << "# Integration weights " << std::endl;
      out << "#----------------------------------------------------" <<
        std::endl;
      for(int ip=0; ip < m_nIntegrationPoints; ip++){  
        out << m_integrationWeights[ip] << "  ";
        out << std::endl;
        sum += m_integrationWeights[ip];
      }

      out << std::endl;
      out << "#Sum = " << sum << std::endl;
      out << "#----------------------------------------------------" <<
        std::endl;
    }
    return;
  }

  /********************************************//**
   * \brief 
   *
   * 
   ***********************************************/
    template<class T_S,class T_P>
    void SurrogatePseudoSpectral<T_S,T_P>::prettyPrintIntegrationWeights( ) const
    {
      if (this->m_comm->rank() == 0)
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
      }
      return;
    }


  /********************************************//**
   * \brief 
   ***********************************************/
    template<class T_S, class T_P> 
      void SurrogatePseudoSpectral<T_S,T_P>::printIndexSet( std::ostream& out ) const
      {
        if (this->m_comm->rank() == 0)
        {
          out << std::endl;
          out << "#====================================================" <<
            std::endl;
          out << "# Index Set" << std::endl;
          out << "#----------------------------------------------------" <<
          std::endl;

          for (unsigned int i=0; i< m_indexSet.size(); i++)
          {
            for (unsigned int j=0; j< m_indexSet[i].size(); j++)
            {
              out << std::setw(5) << m_indexSet[i][j] << " " ;
            }
            out << std::endl;
          }
          out << "#----------------------------------------------------" <<
            std::endl;
        }

      return ;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P> 
    void SurrogatePseudoSpectral<T_S,T_P>::prettyPrintIndexSet( ) const
    {
      if (this->m_comm->rank() == 0)
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

      }
      return ;
    }

  /********************************************//**
   * \brief calculate l2 norm by integration
   ***********************************************/
  template<class T_S, class T_P>
    std::map< std::string, T_P> SurrogatePseudoSpectral<T_S,T_P>::l2Norm(
        std::vector< std::string > solutionNames
        )
    {

      // initialize to zero
      std::map< std::string, T_P> l2norm;
      for (unsigned int n=0; n < solutionNames.size(); n++)
      {
        std::string id = solutionNames[n];
        l2norm.insert( std::pair<std::string,T_P>(
              id, T_P( std::vector<double>(this->m_solSize[id],0.) ) )
            );
      }

      //---------
      //get integration points for higher order quad rule
      std::vector<unsigned int> integrationOrder = this->m_order;
      for(unsigned int i=0; i<integrationOrder.size();i++)
        integrationOrder[i] += 2;

      QuadratureTensorProduct integrationQuadratureRule(
          this->m_parameters, integrationOrder );

      unsigned int dimension = this->m_dimension;
      unsigned int nQuadPoints = integrationQuadratureRule.getNQuadPoints();
      double** quadPoints = integrationQuadratureRule.getQuadPoints();
      double* quadWeights = integrationQuadratureRule.getQuadWeights();



      //----------
      // loop over quad points
      for(unsigned int i=0; i<nQuadPoints; i++)
      {
        T_S paramValues(dimension);
        for(unsigned int j=0; j<paramValues.size(); j++)
          paramValues(j) = quadPoints[i][j];

        std::map<std::string,T_P> pointValue 
          = this->evaluate(solutionNames,paramValues);
        for (unsigned int n=0; n < solutionNames.size(); n++)
        {
          std::string id = solutionNames[n];
          for(unsigned int j=0; j<l2norm[id].size();j++)
            l2norm[id](j) += pointValue[id](j) * pointValue[id](j) 
              * quadWeights[i] ; 
        }

      } // end loop over quad points

      
      // take square root
      for (unsigned int n=0; n < solutionNames.size(); n++)
      {
        std::string id = solutionNames[n];

        for(unsigned int j=0; j<l2norm[id].size();j++)
          l2norm[id](j) =  std::sqrt( l2norm[id](j) );
      }


      return l2norm;
    }


  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
    double SurrogatePseudoSpectral<T_S,T_P>::l2NormDifference(
        SurrogateModel<T_S,T_P>& comparisonModel,
        std::string solutionName )
    {

      // make sure solutionName is present in both models
      // evaluate at integration points
      //
      //
      // initialize to zero
      double l2norm = 0;

      if (this->m_solutionFunction.count(solutionName) == 0) 
      {
        std::cerr << "\n\t ERROR: requested solution name not present in "
          << " base model of l2NormDifference( ... ). \n "
          << std::endl;
        exit(1);
      }
      if (comparisonModel.getSolutionNames().count(solutionName) == 0 )
      {
        std::cerr << "\n\t ERROR: requested solution name not present in "
          << " comparison model of l2NormDifference( ... ). \n "
          << std::endl;
        exit(1);
      }

      //---------
      //get integration points for higher order quad rule
      std::vector<unsigned int> integrationOrder = this->m_order;
      std::vector<unsigned int> comparisonOrder 
        = comparisonModel.getExpansionOrder() ;
      for(unsigned int i=0; i<integrationOrder.size();i++)
        integrationOrder[i] 
          = std::max(integrationOrder[i],comparisonOrder[i]) 
          + 2;

      QuadratureTensorProduct integrationQuadratureRule(
          this->m_parameters, integrationOrder );

      unsigned int dimension = this->m_dimension;
      unsigned int nQuadPoints = integrationQuadratureRule.getNQuadPoints();
      double** quadPoints = integrationQuadratureRule.getQuadPoints();
      double* quadWeights = integrationQuadratureRule.getQuadWeights();



      //----------
      // loop over quad points
      for(unsigned int i=0; i<nQuadPoints; i++)
      {
        T_S paramValues(dimension);
        for(unsigned int j=0; j<paramValues.size(); j++)
          paramValues(j) = quadPoints[i][j];

        T_P diffVec = this->evaluate(solutionName,paramValues);
        diffVec -= comparisonModel.evaluate(solutionName,paramValues);


        l2norm += diffVec.dot( diffVec )  * quadWeights[i] ; 
        /* l2norm += 1  * quadWeights[i] ; */ 

      } // end loop over quad points

      
      // take square root
      l2norm =  std::sqrt( l2norm );


      return l2norm;

    }
  
}
#endif // SURROGATE_PSEUDO_SPECTRAL_H


