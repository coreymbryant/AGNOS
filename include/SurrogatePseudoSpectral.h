
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
  template<class T_S>
    class SurrogatePseudoSpectral : public SurrogateModel<T_S>
  {

    public:

      // single physics function constructors
      SurrogatePseudoSpectral( 
          const Communicator&               comm,
          PhysicsModel<T_S>*                physics,
          const std::vector<Parameter*>     parameters,
          const unsigned int                order 
          );
      SurrogatePseudoSpectral( 
          const Communicator&               comm,
          PhysicsModel<T_S>*                physics,
          const std::vector<Parameter*>     parameters,
          const std::vector<unsigned int>&  order
          );

      virtual ~SurrogatePseudoSpectral( );

      void build( ) ;
      std::map< std::string, NumericVector<Number>*  > computeContribution( 
          PhysicsModel<T_S>*  physics,
          T_S&                integrationPoint, 
          double&             integrationWeight
          );

      using SurrogateModel<T_S>::evaluate; 
      std::map<std::string, std::vector<Number> > evaluate( 
          std::vector< std::string >  solutionNames,
          T_S&                        parameterValues 
          );

      /* using SurrogateModel<T_S>::l2Norm; */
      /* std::map< std::string, NumericVector<Number> > l2Norm( */
      /*   std::vector< std::string > solutionNames */
      /*   ); */
      /* double l2NormDifference( */
      /*     SurrogateModel<T_S>& comparisonModel, */
      /*     std::string solutionName ); */

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
 * \brief Constructor for isotropic order
 ***********************************************/
  template<class T_S>
    SurrogatePseudoSpectral<T_S>::SurrogatePseudoSpectral( 
        const Communicator&               comm,
        PhysicsModel<T_S>*                physics,
        const std::vector<Parameter*> parameters,
        const unsigned int order
        )
      : SurrogateModel<T_S>(comm,physics,parameters,order) 
    {
    }


/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S>
    SurrogatePseudoSpectral<T_S>::SurrogatePseudoSpectral( 
        const Communicator&               comm,
        PhysicsModel<T_S>*                physics,
        const std::vector<Parameter*> parameters,
        const std::vector<unsigned int>& order
        )
      : SurrogateModel<T_S>(comm,physics,parameters,order) 
    {
    }


/********************************************//**
 * \brief Constructor for anisotropic order
 ***********************************************/
  template<class T_S>
    SurrogatePseudoSpectral<T_S>::~SurrogatePseudoSpectral( )
    {
    }

/********************************************//**
 * \brief Get number of integration points
 ***********************************************/
  template<class T_S>
    unsigned int SurrogatePseudoSpectral<T_S>::getNIntegrationPoints( )
    const
  {
    return _nIntegrationPoints;
  }

/********************************************//**
 * \brief Get integration weights
 ***********************************************/
  template<class T_S>
    std::vector<double> SurrogatePseudoSpectral<T_S>::getIntegrationWeights( )
    const
  {
    return _integrationWeights;
  }

/********************************************//**
 * \brief Get the current expansion order
 ***********************************************/
  template<class T_S>
    std::vector<T_S> SurrogatePseudoSpectral<T_S>::getIntegrationPoints( )
    const
  {
    return _integrationPoints;
  }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S> 
    const std::vector< std::vector< unsigned int> > 
    SurrogatePseudoSpectral<T_S>::getIndexSet( ) const
    {
      return _indexSet;
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S> 
    std::vector<double> SurrogatePseudoSpectral<T_S>::evaluateBasis( 
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
  template<class T_S>
    void SurrogatePseudoSpectral<T_S>::build( )
    {
      // This is separated from the routine that actually computes contribution
      // so that we can group surrogate models together later and wll that needs
      // to be defined is build routine based on computeContribution( )
      //
      unsigned int totalNCoeff = this->_totalNCoeff;
      std::vector< std::vector<double> > polyValues;
      polyValues.reserve(totalNCoeff);

      std::set<std::string>::iterator id;
      
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


      // initiaize contribution vector
      std::vector< std::map< std::string, NumericVector<Number>*  > > mySolutionVectors;
      mySolutionVectors.reserve(nPts);

      this->_comm.barrier();


      if(DEBUG)
        std::cout << "order:" << this->_order[0] << std::endl;

      if( this->_comm.rank() == 0)
        std::cout << "     --> Solving at " << _nIntegrationPoints 
          << " integration points " << std::endl;


      // solve for my integration points
      for(unsigned int pt=0; pt < nPts; pt++)
      {
        // compute the contribution
        mySolutionVectors.push_back( computeContribution( 
            this->_physics,
            _integrationPoints[intPtsStart + pt*this->_comm.size() ], 
            _integrationWeights[intPtsStart + pt*this->_comm.size() ]
            ) );

        polyValues.push_back(
            evaluateBasis( 
              this->_indexSet, 
              _integrationPoints[intPtsStart + pt*this->_comm.size()]
              ) 
            );
      }

      this->_comm.barrier();

      // need to know size of solution vector on all processes (in case some
      // don't have any work but will still be called in reduce operation)
      /* unsigned int tempSize; */
      /* this->_solSize.clear(); */
      /* for (id=this->_solutionFunction.begin(); */
      /*     id!=this->_solutionFunction.end(); id++) */
      /* { */
      /*   if (this->_comm.rank() == 0) */
      /*     tempSize = mySolutionVectors[0][id->first].size(); */

      /*   this->_comm.broadcast(tempSize); */

      /*   this->_solSize.insert( std::pair<std::string,unsigned int>( */
      /*         id->first,tempSize ) ); */

      /* } */

      // WAIT FOR ALL PROCESSES TO CATCH UP
      this->_comm.barrier();

      if( this->_comm.rank() == 0)
        std::cout << "     --> Computing " << totalNCoeff << " coefficients" <<
          std::endl;


      // --------------
      // sum my solutions weighted by polyValue
      // localize so that we can send it to other groups
      // sum all processes contributions
      // save if its one of my coefficients
      // --------------
      
      // clear old coefficients
      this->_coefficients.clear();

      //-- loop through sols
      for (id=this->_solutionNames.begin();
          id!=this->_solutionNames.end(); ++id)
      {
        // temporary storage for my coefficients for this sol
        std::vector< std::vector<Number> > solCoefficientVectors;
        solCoefficientVectors.reserve(nCoeffs);

        for (unsigned int c=0; c < totalNCoeff; c++)
        {
          std::vector<Number> sendContrib;

          NumericVector<Number>* myContrib =
            mySolutionVectors[0][*id]->zero_clone().release();

          // compute current rank's contribution
          for(unsigned int pt=0; pt < nPts; pt++)
          {
            NumericVector<Number>* temp =
              mySolutionVectors[pt][*id]->clone().release();
              temp->scale(polyValues[pt][c]);
            myContrib->add( *temp );
          }

          myContrib->localize( sendContrib );

          // sum all contributions
          this->_comm.sum( sendContrib );

          // store if its one of mine
          if ( (c % this->_comm.size()) == this->_comm.rank() )
            solCoefficientVectors.push_back( sendContrib ) ;

        } // c

        if (solCoefficientVectors.size() > 0)
        {
          this->_coefficients.insert(
              std::pair< std::string, std::vector< std::vector<Number> > >(*id,
                solCoefficientVectors) ) ;
        }

      } // id

      return;
    } 

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S>
    std::map< std::string, NumericVector<Number>*  >
    SurrogatePseudoSpectral<T_S>::computeContribution(
        PhysicsModel<T_S>*  physics,
        T_S&                integrationPoint, 
        double&             integrationWeight
        )
    {

      if(DEBUG)
        std::cout << "test: computeContribution( ) beginning" << std::endl;
      

      std::map< std::string, NumericVector<Number>*  > contrib;

      
      if(DEBUG)
        std::cout << "test: computeContribution( ) pre compute()" << std::endl;

      // get solution for this integration point
      physics->compute( integrationPoint, contrib );

      if(DEBUG)
        std::cout << "test: computeContribution( ) post compute()" << std::endl;

      std::set<std::string>::iterator id = this->_solutionNames.begin();
      for(;id!=this->_solutionNames.end();++id)
        contrib[*id]->scale( integrationWeight );


      if(DEBUG)
        std::cout << "test: computeContribution( ) ending" << std::endl;

      return contrib;
    }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  template<class T_S>
    std::map<std::string, std::vector<Number> > SurrogatePseudoSpectral<T_S>::evaluate( 
        std::vector< std::string > solutionNames,
        T_S& parameterValues /**< parameter values to evaluate*/
        )
    {
      unsigned int totalNCoeff = this->_totalNCoeff;
      std::map< std::string, std::vector<Number> > surrogateValue;

      // get polyValues
      std::vector<double> polyValues =
        evaluateBasis(this->_indexSet,parameterValues) ;

      // get starting coeff for current rank
      unsigned int coeffStart = std::min(this->_comm.rank(), totalNCoeff-1);

      // loop over all solutions requested
      for (unsigned int n=0; n < solutionNames.size(); n++)
      {
        std::string id = solutionNames[n];

        std::vector<Number> sendContrib(this->_coefficients[id][0].size(),0.);

        for(unsigned int c=0; c < this->_coefficients[id].size(); c++)
        {
          // compute current rank's contribution
          for(unsigned int i=0; i < this->_coefficients[id][c].size(); i++)
            sendContrib[i] += this->_coefficients[id][c][i] *
              polyValues[coeffStart + c*(this->_comm.size())];

        } // coefficients

        // sum all contributions
        this->_comm.sum( sendContrib );

        surrogateValue.insert(
            std::pair< std::string, std::vector<Number>  >(id, sendContrib ) 
            ) ;

      } // solNames

      return surrogateValue;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template< class T_S>
  void SurrogatePseudoSpectral<T_S>::printIntegrationPoints( 
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
  template< class T_S>
  void SurrogatePseudoSpectral<T_S>::prettyPrintIntegrationPoints( ) const
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
  template< class T_S>
  void SurrogatePseudoSpectral<T_S>::printIntegrationWeights( 
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
    template< class T_S>
    void SurrogatePseudoSpectral<T_S>::prettyPrintIntegrationWeights( ) const
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
    template<class T_S> 
      void SurrogatePseudoSpectral<T_S>::printIndexSet( std::ostream& out ) const
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
  template<class T_S> 
    void SurrogatePseudoSpectral<T_S>::prettyPrintIndexSet( ) const
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
  /* template<class T_S> */
  /*   std::map< std::string, NumericVector<Number> > SurrogatePseudoSpectral<T_S>::l2Norm( */
  /*       std::vector< std::string > solutionNames */
  /*       ) */
  /*   { */

  /*     // initialize to zero */
  /*     std::map< std::string, NumericVector<Number> > l2norm; */
  /*     for (unsigned int n=0; n < solutionNames.size(); n++) */
  /*     { */
  /*       std::string id = solutionNames[n]; */
  /*       l2norm.insert( std::pair<std::string,NumericVector<Number> >( */
  /*             id, NumericVector<Number> ( std::vector<double>(this->_solSize[id],0.) ) ) */
  /*           ); */
  /*     } */

  /*     //--------- */
  /*     //get integration points for higher order quad rule */
  /*     std::vector<unsigned int> integrationOrder = this->_order; */
  /*     for(unsigned int i=0; i<integrationOrder.size();i++) */
  /*       integrationOrder[i] += 2; */

  /*     QuadratureTensorProduct integrationQuadratureRule( */
  /*         this->_parameters, integrationOrder ); */

  /*     unsigned int dimension = this->_dimension; */
  /*     unsigned int nQuadPoints = integrationQuadratureRule.getNQuadPoints(); */
  /*     double** quadPoints = integrationQuadratureRule.getQuadPoints(); */
  /*     double* quadWeights = integrationQuadratureRule.getQuadWeights(); */



  /*     //---------- */
  /*     // loop over quad points */
  /*     for(unsigned int i=0; i<nQuadPoints; i++) */
  /*     { */
  /*       T_S paramValues(dimension); */
  /*       for(unsigned int j=0; j<paramValues.size(); j++) */
  /*         paramValues(j) = quadPoints[i][j]; */

  /*       std::map<std::string,NumericVector<Number> > pointValue */ 
  /*         = this->evaluate(solutionNames,paramValues); */
  /*       for (unsigned int n=0; n < solutionNames.size(); n++) */
  /*       { */
  /*         std::string id = solutionNames[n]; */
  /*         for(unsigned int j=0; j<l2norm[id].size();j++) */
  /*           l2norm[id](j) += pointValue[id](j) * pointValue[id](j) */ 
  /*             * quadWeights[i] ; */ 
  /*       } */

  /*     } // end loop over quad points */

      
  /*     // take square root */
  /*     for (unsigned int n=0; n < solutionNames.size(); n++) */
  /*     { */
  /*       std::string id = solutionNames[n]; */

  /*       for(unsigned int j=0; j<l2norm[id].size();j++) */
  /*         l2norm[id](j) =  std::sqrt( l2norm[id](j) ); */
  /*     } */


  /*     return l2norm; */
  /*   } */


  /********************************************//**
   * \brief 
   ***********************************************/
  /* template< class T_S> */
  /*   double SurrogatePseudoSpectral<T_S>::l2NormDifference( */
  /*       SurrogateModel<T_S>& comparisonModel, */
  /*       std::string solutionName ) */
  /*   { */

  /*     // make sure solutionName is present in both models */
  /*     // evaluate at integration points */
  /*     // */
  /*     // */
  /*     // initialize to zero */
  /*     double l2norm = 0; */

  /*     if (this->_solutionFunction.count(solutionName) == 0) */ 
  /*     { */
  /*       std::cerr << "\n\t ERROR: requested solution name not present in " */
  /*         << " base model of l2NormDifference( ... ). \n " */
  /*         << std::endl; */
  /*       exit(1); */
  /*     } */
  /*     if (comparisonModel.getSolutionNames().count(solutionName) == 0 ) */
  /*     { */
  /*       std::cerr << "\n\t ERROR: requested solution name not present in " */
  /*         << " comparison model of l2NormDifference( ... ). \n " */
  /*         << std::endl; */
  /*       exit(1); */
  /*     } */

  /*     //--------- */
  /*     //get integration points for higher order quad rule */
  /*     std::vector<unsigned int> integrationOrder = this->_order; */
  /*     std::vector<unsigned int> comparisonOrder */ 
  /*       = comparisonModel.getExpansionOrder() ; */
  /*     for(unsigned int i=0; i<integrationOrder.size();i++) */
  /*       integrationOrder[i] */ 
  /*         = std::max(integrationOrder[i],comparisonOrder[i]) */ 
  /*         + 2; */

  /*     QuadratureTensorProduct integrationQuadratureRule( */
  /*         this->_parameters, integrationOrder ); */

  /*     unsigned int dimension = this->_dimension; */
  /*     unsigned int nQuadPoints = integrationQuadratureRule.getNQuadPoints(); */
  /*     double** quadPoints = integrationQuadratureRule.getQuadPoints(); */
  /*     double* quadWeights = integrationQuadratureRule.getQuadWeights(); */



  /*     //---------- */
  /*     // loop over quad points */
  /*     for(unsigned int i=0; i<nQuadPoints; i++) */
  /*     { */
  /*       T_S paramValues(dimension); */
  /*       for(unsigned int j=0; j<paramValues.size(); j++) */
  /*         paramValues(j) = quadPoints[i][j]; */

  /*       NumericVector<Number>  diffVec = this->evaluate(solutionName,paramValues); */
  /*       diffVec -= comparisonModel.evaluate(solutionName,paramValues); */


  /*       l2norm += diffVec.dot( diffVec )  * quadWeights[i] ; */ 
  /*       /1* l2norm += 1  * quadWeights[i] ; *1/ */ 

  /*     } // end loop over quad points */

      
  /*     // take square root */
  /*     l2norm =  std::sqrt( l2norm ); */


  /*     return l2norm; */

  /*   } */
  
}
#endif // SURROGATE_PSEUDO_SPECTRAL_H


