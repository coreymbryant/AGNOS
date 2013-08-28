
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

      /** Constructor */
      SurrogatePseudoSpectral( 
          const Communicator&               comm,
          PhysicsModel<T_S,T_P>*                physics,
          const std::vector<Parameter*>     parameters,
          const std::vector<unsigned int>&  order
          );


      /** Secondary Constructor. 
       *  ***DO NOT MISTAKE FOR A COPY CONSTRUCTOR ***
       *  Intended use is for constructing a secondary surrogate model using the
       *  primary model as an evaluating object in the build routine. 
       *  If additional inputs are defined it will
       * construct a new surrogate increasing the order and using
       * primarySurrogate to perform evaluations in the constructions */
      SurrogatePseudoSpectral( 
          const SurrogateModel<T_S,T_P>* primarySurrogate, 
          unsigned int increaseOrder = 0,
          unsigned int multiplyOrder = 1,
          std::set<std::string> evaluateSolutions = std::set<std::string>(),
          std::set<std::string> computeSolutions = std::set<std::string>()
          );

      virtual ~SurrogatePseudoSpectral( );

      void build( ) ;
      std::map< std::string, T_P > computeContribution( 
          PhysicsModel<T_S,T_P>*   physics,
          unsigned int index
          );

      using SurrogateModel<T_S,T_P>::evaluate; 
      std::map<std::string, T_P> evaluate( 
          std::set< std::string >  solutionNames,
          T_S&                        parameterValues 
          ) const ;

      using SurrogateModel<T_S,T_P>::l2Norm;
      std::map< std::string, T_P> l2Norm(
        std::set< std::string > solutionNames
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

      unsigned int              _nIntegrationPoints;
      std::vector<T_S>          _integrationPoints ;
      std::vector<double>       _integrationWeights ;

      std::vector< std::vector<unsigned int> > _indexSet;



  };



/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        const Communicator&               comm,
        PhysicsModel<T_S,T_P>*                physics,
        const std::vector<Parameter*>     parameters,
        const std::vector<unsigned int>& order
        )
      : SurrogateModel<T_S,T_P>( comm, physics, parameters, order) 
    {
    }


/********************************************//*
 * \brief Secondary constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogatePseudoSpectral<T_S,T_P>::SurrogatePseudoSpectral( 
        const SurrogateModel<T_S,T_P>* primarySurrogate, 
        unsigned int increaseOrder ,
        unsigned int multiplyOrder ,
        std::set<std::string> evaluateSolutions,
        std::set<std::string> computeSolutions
        )
      : SurrogateModel<T_S,T_P>(primarySurrogate, increaseOrder, multiplyOrder,
          evaluateSolutions, computeSolutions)
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
    return _nIntegrationPoints;
  }

/********************************************//**
 * \brief Get integration weights
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<double> SurrogatePseudoSpectral<T_S,T_P>::getIntegrationWeights( )
    const
  {
    return _integrationWeights;
  }

/********************************************//**
 * \brief Get the current expansion order
 ***********************************************/
  template<class T_S, class T_P>
    std::vector<T_S> SurrogatePseudoSpectral<T_S,T_P>::getIntegrationPoints( )
    const
  {
    return _integrationPoints;
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
      return _indexSet;
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
        for(unsigned int dir=0; dir < this->_dimension; dir++)
        {
          basisValues[id] 
            *= this->_parameters[dir]->evalBasisPoly( 
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
      unsigned int totalNCoeff = this->_totalNCoeff;
      std::vector< std::vector<double> > polyValues;
      polyValues.reserve(totalNCoeff);
      std::set< std::string >::iterator id;
      
      // output progress
      if( this->_comm.rank() == 0)
      {
        std::cout << std::endl;
        std::cout << " Starting build of surrogate model for: { " ;
        for (id=this->_solutionNames.begin();
            id!=this->_solutionNames.end(); id++)
        {
          if (id!=this->_solutionNames.begin())
            std::cout << " , " ;
          std::cout << *id;
        }
        std::cout << " }" << std::endl;
      }
      
      // ------------------------------------
      // Int points and coeff are currently separated in case we want to
      // implememnt higher order integration rule in the future
      // assign the appropriate int points to each processor
      unsigned int nPts  
        =  _nIntegrationPoints / this->_comm.size()  
        + ( this->_comm.rank() < (_nIntegrationPoints % this->_comm.size() ) ) ;
      unsigned int intPtsStart = std::min(this->_comm.rank(), _nIntegrationPoints-1);



      // assign the appropriate coeff to each processor
      unsigned int nCoeffs  
        =  totalNCoeff / this->_comm.size()  
        + ( this->_comm.rank() < (totalNCoeff % this->_comm.size() ) ) ;
      unsigned int coeffStart = std::min(this->_comm.rank(), totalNCoeff-1);
      // ------------------------------------
      
      
      // ------------------------------------
      // Get primalSurrogate evaluations if needed
      this->_primalEvaluations.clear();
      if ( ! this->_evalNames.empty() )
      {
        std::cout << "_integrationPoints.size():"  << _integrationPoints.size() << std::endl;
        for(unsigned int pt=0; pt < _nIntegrationPoints; pt++)
        {
          std::cout << "pt: " << pt << std::endl;
          T_S integrationPoint = _integrationPoints[pt];
          this->_primalEvaluations.push_back( 
              this->_evalSurrogate->evaluate( this->_evalNames, integrationPoint  )
              );
        }
      }
      
      // ------------------------------------





      this->_comm.barrier();


      // ------------------------------------
      // Solve for solution at integration points
      if( this->_comm.rank() == 0)
        std::cout << "     --> Solving at " << _nIntegrationPoints 
          << " integration points " << std::endl;
      
      // initiaize contribution vector
      std::vector< std::map< std::string, T_P > > myContribs;
      myContribs.reserve(nPts);


      // solve for my integration points
      for(unsigned int pt=0; pt < nPts; pt++)
      {
        if (AGNOS_DEBUG)
          std::cout << "test: beginning of pt" << std::endl;



        if (AGNOS_DEBUG)
          std::cout << "test: pt (pre computeContribution)" << std::endl;

        // compute the contribution
        myContribs.push_back( computeContribution( 
            this->_physics,
            intPtsStart + pt*this->_comm.size()
            ) );

        if (AGNOS_DEBUG)
          std::cout << "test: pt (post computeContribution)" << std::endl;

        polyValues.push_back(
            evaluateBasis( this->_indexSet, _integrationPoints[intPtsStart +
              pt*this->_comm.size()]) 
            );

        if (AGNOS_DEBUG)
          std::cout << "test: end of pt" << std::endl;
      }

        std::cout << "rank: " << this->_comm.rank() << std::endl;
        std::cout << "size: " << this->_comm.size() << std::endl;
      this->_comm.barrier();


      // need to know size of solution vector on all processes (in case some
      // don't have any work but will still be called in reduce operation)
      unsigned int tempSize;
      this->_solSize.clear();
      for (id=this->_solutionNames.begin();
          id!=this->_solutionNames.end(); id++)
      {
        if (this->_comm.rank() == 0)
          tempSize = myContribs[0][*id].size();

        if(AGNOS_DEBUG)
          std::cout << "surrogate build -- solSize[" << *id << "]:" <<
            tempSize << std::endl;

        this->_comm.broadcast(tempSize);

        if(AGNOS_DEBUG)
          std::cout << "surrogate build --(post_broadcast) solSize[" << *id << "]:" <<
            tempSize << std::endl;

        this->_solSize.insert( std::pair<std::string,unsigned int>(
              *id,tempSize ) );

      }



      // WAIT FOR ALL PROCESSES TO CATCH UP
      this->_comm.barrier();
      if( this->_comm.rank() == 0)
        std::cout << "     --> Computing " << totalNCoeff << " coefficients" <<
          std::endl;


      // --------------
      // compute my contribution to each coefficient
      // sum all processes contributions
      // save if its one of my coefficients
      // --------------
      
      // clear old coefficients
      this->_coefficients.clear();

      if(AGNOS_DEBUG)
        std::cout << "surrogate build before computing coeffs" << std::endl ;

      //-- loop through sols
      for (id=this->_solutionNames.begin();
          id!=this->_solutionNames.end(); id++)
      {
        if(AGNOS_DEBUG)
          std::cout << "surrogate build computing coeffs: " << *id << std::endl ;
        // temporary storage for my coefficients for this sol
        std::vector<T_P> solCoefficientVectors;
        solCoefficientVectors.reserve(nCoeffs);

        for (unsigned int c=0; c < totalNCoeff; c++)
        {
          std::vector<double> sumContrib(this->_solSize[*id],0.);
          
          // compute current rank's contribution
          for(unsigned int pt=0; pt < nPts; pt++)
            for(unsigned int i=0; i <  sumContrib.size(); i++)
              sumContrib[i] += myContribs[pt][*id](i) *
                polyValues[pt][c];

          // sum all contributions
          this->_comm.sum( sumContrib );

          // store if its one of mine
          if ( (c % this->_comm.size()) == this->_comm.rank() )
            solCoefficientVectors.push_back( T_P(sumContrib) ) ;

        } // c

        if (solCoefficientVectors.size() > 0)
        {
          this->_coefficients.insert(
              std::pair< std::string, std::vector<T_P> >(
                *id,
                solCoefficientVectors
                ) 
              );
        }

      } // id

        if(AGNOS_DEBUG)
          std::cout << "leaving surrogate build " << std::endl ;
      return;
    } 

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    std::map< std::string, T_P >
    SurrogatePseudoSpectral<T_S,T_P>::computeContribution(
        PhysicsModel<T_S,T_P>*  physics,
        unsigned int index
        )
    {

      if (AGNOS_DEBUG)
        std::cout << "test: computeContribution( ) beginning" << std::endl;
      

      // get data for this integration point
      std::map< std::string, T_P > contrib ;
      if ( ! this->_primalEvaluations.empty() )
        contrib = this->_primalEvaluations[index];
      T_S integrationPoint = _integrationPoints[index];
      double integrationWeight = _integrationWeights[index];

      
      if (AGNOS_DEBUG)
        std::cout << "test: computeContribution( ) pre compute()" << std::endl;

      // get solution for this integration point
      physics->compute( integrationPoint, contrib );


      // clear out any evaluations that were provided by evalSurrogate, and that
      // aren't in solutionNames
      std::set<std::string>::iterator evalName = this->_evalNames.begin() ;
      for( ; evalName != this->_evalNames.end(); evalName++)
      {
        if ( !this->_solutionNames.count(*evalName) ) 
          contrib.erase( *evalName );
      }
      

      if (AGNOS_DEBUG)
      {
        std::cout << "test: computeContribution( ) post compute()" << std::endl;
      }

      std::set<std::string>::iterator id = this->_solutionNames.begin();
      for(;id!=this->_solutionNames.end();++id)
        contrib[*id].scale( integrationWeight );


      if (AGNOS_DEBUG)
        std::cout << "test: computeContribution( ) ending" << std::endl;
      return contrib;
    }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  template<class T_S, class T_P>
    std::map<std::string, T_P> SurrogatePseudoSpectral<T_S,T_P>::evaluate( 
        std::set< std::string > solutionNames,
        T_S& parameterValues /**< parameter values to evaluate*/
        ) const
    {
      if(AGNOS_DEBUG)
        std::cout << "entering surrogate evaluate" << std::endl;
      unsigned int totalNCoeff = this->_totalNCoeff;
      std::map< std::string, T_P> surrogateValue;

      // get polyValues
      std::vector<double> polyValues =
        evaluateBasis(this->_indexSet,parameterValues) ;

      // get starting coeff for current rank
      unsigned int coeffStart = std::min(this->_comm.rank(), totalNCoeff-1);

      // get references to system variables. We can't iterate directly since
      // this is a const member funciton
      std::map<std::string,unsigned int> solSize = this->getSolSize();
      std::map<std::string,std::vector<T_P> > coefficients =
        this->getLocalCoefficients();

      // loop over all solutions requested
      std::set<std::string>::iterator id = solutionNames.begin();
      for (; id != solutionNames.end(); id++)
      {

        std::vector<double> sumContrib(solSize[*id],0.);
        if(AGNOS_DEBUG)
          std::cout << "evaluating surrogate model for: " << *id <<
            "with size:" << sumContrib.size() << std::endl;

        for(unsigned int c=0; c < coefficients[*id].size(); c++)
        {
          // compute current rank's contribution
          for(unsigned int i=0; i < coefficients[*id][c].size(); i++)
            sumContrib[i] += coefficients[*id][c](i) *
              polyValues[coeffStart + c*(this->_comm.size())];

        } // coefficients

        // sum all contributions
        this->_comm.sum( sumContrib );

        surrogateValue.insert(
            std::pair< std::string, T_P >(*id, T_P(sumContrib) ) 
            ) ;

      } // solNames

      if(AGNOS_DEBUG)
        std::cout << "leaving surrogate evaluate" << std::endl;
      return surrogateValue;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
  void SurrogatePseudoSpectral<T_S,T_P>::printIntegrationPoints( 
      std::ostream& out ) const
  {
    if (this->_comm.rank() == 0)
    {
      out << std::endl;
      out << "#====================================================" <<
        std::endl;
      out << "# Integration points " << std::endl;
      out << "#----------------------------------------------------" <<
        std::endl;
      for(int ix=0; ix < _nIntegrationPoints ; ix++)  
      {
        for(int iy=0; iy < this->_dimension; iy++)
        {
          out << std::scientific << std::setprecision(5) << std::setw(12)
            << _integrationPoints[ix](iy) << "  ";
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
    if (this->_comm.rank() == 0)
    {
      std::cout << std::endl;
      std::cout << "====================================================" <<
        std::endl;
      std::cout << " Integration points " << std::endl;
      std::cout << "----------------------------------------------------" <<
        std::endl;
      std::cout << "   \\ x        " << std::endl;
      std::cout << "    \\" ;
      for(unsigned int dim=0; dim < this->_dimension; dim++)
        std::cout << std::setw(12) << "x_" << dim << " " ;
      std::cout << std::endl;
      std::cout << "  id \\  " << std::endl;
      std::cout << "----------------------------------------------------" <<
        std::endl;
      for(int ix=0; ix < _nIntegrationPoints ; ix++)  
      {
        std::cout << std::setw(5) << ix << " |   ";
        for(int iy=0; iy < this->_dimension; iy++)
        {
          std::cout << std::scientific << std::setprecision(5) << std::setw(12)
            << _integrationPoints[ix](iy) << "  ";
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
    if (this->_comm.rank() == 0)
    {
      double sum = 0.0;
      out << std::endl;
      out << "#====================================================" <<
        std::endl;
      out << "# Integration weights " << std::endl;
      out << "#----------------------------------------------------" <<
        std::endl;
      for(int ip=0; ip < _nIntegrationPoints; ip++){  
        out << _integrationWeights[ip] << "  ";
        out << std::endl;
        sum += _integrationWeights[ip];
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
      if (this->_comm.rank() == 0)
      {
        double sum = 0.0;
        std::cout << std::endl;
        std::cout << "====================================================" <<
          std::endl;
        std::cout << " Integration weights " << std::endl;
        std::cout << "----------------------------------------------------" <<
          std::endl;
        for(int ip=0; ip < _nIntegrationPoints; ip++){  
          std::cout << std::setw(5) << ip << "   |   ";
          std::cout << _integrationWeights[ip] << "  ";
          std::cout << std::endl;
          sum += _integrationWeights[ip];
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
        if (this->_comm.rank() == 0)
        {
          out << std::endl;
          out << "#====================================================" <<
            std::endl;
          out << "# Index Set" << std::endl;
          out << "#----------------------------------------------------" <<
          std::endl;

          for (unsigned int i=0; i< _indexSet.size(); i++)
          {
            for (unsigned int j=0; j< _indexSet[i].size(); j++)
            {
              out << std::setw(5) << _indexSet[i][j] << " " ;
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
      if (this->_comm.rank() == 0)
      {
        std::cout << std::endl;
        std::cout << "====================================================" <<
          std::endl;
        std::cout << " Index Set" << std::endl;
        std::cout << "----------------------------------------------------" <<
        std::endl;
      std::cout << "   \\ dir      " << std::endl;
      std::cout << "    \\        " ;
      for(unsigned int dim=0; dim < this->_dimension; dim++)
        std::cout << std::setw(4) << "xi_" << dim << " " ;
      std::cout << std::endl;
      std::cout << "  id \\  " << std::endl;
      std::cout << "----------------------------------------------------" <<
        std::endl;
        for (unsigned int i=0; i< _indexSet.size(); i++)
        {
        std::cout << std::setw(5) << i << "   |   ";
          for (unsigned int j=0; j< _indexSet[i].size(); j++)
          {
            std::cout << std::setw(5) << _indexSet[i][j] << " " ;
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
        std::set< std::string > solutionNames
        )
    {

      // initialize to zero
      std::map< std::string, T_P> l2norm;
      std::set<std::string>::iterator id = solutionNames.begin();
      for (; id != solutionNames.end(); id++)
      {
        l2norm.insert( std::pair<std::string,T_P>(
              *id, T_P( std::vector<double>(this->_solSize[*id],0.) ) )
            );
      }

      //---------
      //get integration points for higher order quad rule
      std::vector<unsigned int> integrationOrder = this->_order;
      for(unsigned int i=0; i<integrationOrder.size();i++)
        integrationOrder[i] += 2;

      QuadratureTensorProduct integrationQuadratureRule(
          this->_parameters, integrationOrder );

      unsigned int dimension = this->_dimension;
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
        for (id=solutionNames.begin(); id != solutionNames.end(); id++)
        {
          for(unsigned int j=0; j<l2norm[*id].size();j++)
            l2norm[*id](j) += pointValue[*id](j) * pointValue[*id](j) 
              * quadWeights[i] ; 
        }

      } // end loop over quad points

      
      // take square root
      for (id=solutionNames.begin(); id != solutionNames.end(); id++)
      {

        for(unsigned int j=0; j<l2norm[*id].size();j++)
          l2norm[*id](j) =  std::sqrt( l2norm[*id](j) );
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

      if (this->_solutionNames.count(solutionName) == 0) 
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
      std::vector<unsigned int> integrationOrder = this->_order;
      std::vector<unsigned int> comparisonOrder 
        = comparisonModel.getExpansionOrder() ;
      for(unsigned int i=0; i<integrationOrder.size();i++)
        integrationOrder[i] 
          = std::max(integrationOrder[i],comparisonOrder[i]) 
          + 2;

      QuadratureTensorProduct integrationQuadratureRule(
          this->_parameters, integrationOrder );

      unsigned int dimension = this->_dimension;
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


