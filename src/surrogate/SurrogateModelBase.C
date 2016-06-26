
#include "SurrogateModelBase.h"
#include "QuadratureTensorProduct.h"
#include <gsl/gsl_rng.h>
namespace AGNOS
{

  /********************************************//**
   * \brief Constructor
   ***********************************************/
  template<class T_S,class T_P>
    SurrogateModelBase<T_S,T_P>::SurrogateModelBase(
        const Communicator&               comm,
        const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
        const std::vector<unsigned int>&          order,
        std::set<std::string> computeSolutions 
        )
    :
      _comm(comm),
      _parameters(parameters),
      _dimension( parameters.size() ),
      _order(order),
      _solutionNames(computeSolutions),
      _physics(NULL),
      _physicsGroup(0),
      _nPhysicsGroups(1),
      _groupRank(0) 
    {

      if (_order.size() != parameters.size() )
      {
        std::cout 
          << std::endl
          << "\tERROR:"
          << " order vector dimension does not match number of parameters"
          << std::endl
          << std::endl;
        assert(0);
      }
      // we have to set the order for constant parameters and make sure they are
      // set to 0th order
      else
        for( unsigned int i=0; i<_parameters.size(); i++ )
          if ( _parameters[i]->type() == CONSTANT )
            if (_order[i] != 0)
            {
              std::cout 
                << "WARNING: forcing CONSTANT parameter (" 
                << i 
                << ") order to 0th order \n" ;
              _order[i] = 0;
            }
    }


/********************************************//**
 * \brief Default destructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModelBase<T_S,T_P>::~SurrogateModelBase( )
    {
      _coefficients.clear();
    }

/********************************************//**
 * \brief set parameters
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModelBase<T_S,T_P>::setParameters( 
        std::vector<std::shared_ptr<AGNOS::Parameter> >& parameters)
    {
      _parameters = parameters ;
      _dimension = parameters.size();
      return;
    }

/********************************************//**
 * \brief return all coefficients to rank 0, to save on communication we do not
 * broadcast to all processes
 ***********************************************/
  template<class T_S, class T_P> 
    const std::map< std::string, LocalMatrix >
    SurrogateModelBase<T_S,T_P>::getCoefficients( ) const
    {
      if(AGNOS_DEBUG)
        std::cout << "DEBUG: entering getCoefficients routine" << std::endl;
     
      std::map< std::string, LocalMatrix > allCoefficients;
      unsigned int myRank = this->_physicsGroup;

      if (this->_groupRank==0)
      {


        // loop over all solution names
        std::set<std::string>::iterator id = _solutionNames.begin();
        for (; id!=_solutionNames.end(); ++id)
        {
          if(AGNOS_DEBUG)
            std::cout << "DEBUG: getCoefficients routine: id: " << *id << std::endl;

          unsigned int solSize 
            = const_cast<SurrogateModelBase<T_S,T_P>*>(this)->_solSize[*id];
          LocalMatrix solCoefficients(this->_totalNCoeff,solSize);

          if (myRank == 0)
          {
            /* std::cout << "rank 0" << std::endl; */
            /* std::cout << "rank " << myRank << " _coeffIndices.front(): " << */
            /*   this->_coeffIndices.front() << std::endl; */
            /* std::cout << "rank " << myRank << " _coeffIndices.back(): " << */
            /*   this->_coeffIndices.back() << std::endl; */
            /* std::cout << "rank " << myRank << " _coefficients.row_start(): " << */
            /*   this->_coefficients[*id]->row_start() << std::endl; */
            /* std::cout << "rank " << myRank << " _coefficients.row_stop(): " << */
            /*   this->_coefficients[*id]->row_stop() << std::endl; */

            // get my coeffs
            for (unsigned int i=this->_coeffIndices.front(); 
                  i<this->_coeffIndices.back()+1; i++)
              for (unsigned int j=0; j<solSize; j++)
              {
                /* std::cout<<"solCoeff("<<i<<","<<j<<"):" */
                /*   << (*this->_coefficients[*id])(i,j) << std::endl; */
                solCoefficients(i,j) 
                  = (*const_cast<SurrogateModelBase<T_S,T_P>*>(this)
                      ->_coefficients[*id])(i,j)  ;
              }


            //receive other coeffs
            for (unsigned int p=1; p<this->_comm.size(); p++)
            {
              /* std::cout << "rank 0: recv from rank: " << p << std::endl; */
              std::vector<unsigned int> jbuf;
              std::vector<double> cbuf;

              this->_comm.receive(p, jbuf);
              this->_comm.receive(p, cbuf);

              /* std::cout << "rank 0: recv from rank: " << p << " jbuf.size: " << */
              /*   jbuf.size()<< std::endl; */
              /* std::cout << "rank 0: recv from rank: " << p << " cbuf.size: " << */
              /*   cbuf.size()<< std::endl; */

              if (!jbuf.empty())
              {
                for (unsigned int i=jbuf.front(); i<jbuf.back()+1; i++)
                  for (unsigned int j=0; j<solSize; j++)
                  {
                    unsigned int index = 
                      (i-jbuf.front()) * solSize  // due to j index
                      + j ;
                    /* std::cout<<"solCoeff("<<i<<","<<j<<"):" */
                    /*   << cbuf[index] << std::endl; */
                    solCoefficients(i,j) = cbuf[index]  ;
                  } // loop over sol index
              } // if jbuf not empty
            } // for each proc
          } // if rank 0
          else // send my coeffs to proc 0
          {
            /* std::cout << "rank " << myRank << std::endl; */
            std::vector<unsigned int> jbuf;
            std::vector<double> cbuf;

            /* std::cout << "rank " << myRank << " _coeffIndices.size: " << */
            /*   this->_coeffIndices.size() << std::endl; */

            if ( !this->_coeffIndices.empty() )
            {
            /* std::cout << "rank " << myRank << " _coeffIndices.front(): " << */
            /*   this->_coeffIndices.front() << std::endl; */
            /* std::cout << "rank " << myRank << " _coeffIndices.back(): " << */
            /*   this->_coeffIndices.back() << std::endl; */
            /* std::cout << "rank " << myRank << " _coefficients.m(): " << */
            /*   this->_coefficients[*id]->m() << std::endl; */
            /* std::cout << "rank " << myRank << " _coefficients.n(): " << */
            /*   this->_coefficients[*id]->n() << std::endl; */
            /* std::cout << "rank " << myRank << " _coefficients.row_start(): " << */
            /*   this->_coefficients[*id]->row_start() << std::endl; */
            /* std::cout << "rank " << myRank << " _coefficients.row_stop(): " << */
            /*   this->_coefficients[*id]->row_stop() << std::endl; */
              for (unsigned int i=this->_coeffIndices.front(); 
                    i<this->_coeffIndices.back()+1; i++)
              {
                jbuf.push_back(i) ;

                for (unsigned int j=0; j<solSize; j++)
                {
                  /* std::cout << "rank " << myRank << " _coefficients(i,j): " << */
                  /*   (*this->_coefficients[*id])(i,j) << std::endl; */
                  cbuf.push_back( 
                      (*const_cast<SurrogateModelBase<T_S,T_P>*>(this)
                       ->_coefficients[*id])(i,j) 
                      )  ;
                }

              } // for j in index set
            }

            /* std::cout << "rank " << myRank << " jbuf.size: " << jbuf.size() << std::endl; */
            /* std::cout << "rank " << myRank << " cbuf.size: " << cbuf.size() << std::endl; */

            this->_comm.send(0,jbuf);
            this->_comm.send(0,cbuf);


          } // if not proc 0


          allCoefficients.insert(
              std::pair<std::string,LocalMatrix>(*id,solCoefficients) 
              ) ;

          } // end id
      }

      if(AGNOS_DEBUG)
        std::cout << "DEBUG: leaving getCoefficients routine" << std::endl;
      return allCoefficients;
    }

/********************************************//**
 * \brief print coefficient values
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModelBase<T_S,T_P>::printCoefficients( 
        std::vector<std::string> solutionNames,
        std::ostream& out ) 
    {

      if ( this->_groupRank == 0)
      {
        std::map<std::string,LocalMatrix > coefficients =
          this->getCoefficients();

        if ( this->_physicsGroup == 0)
        {

          for (unsigned int i=0; i < solutionNames.size(); i++)
          {
            std::string id = solutionNames[i];
            out << "#" << std::string(75,'-') << std::endl;
            out << "#" << "\t Solution: " << id << std::endl;
            out << "#" << std::string(75,'-') << std::endl;

            for(unsigned int i=0; i<coefficients[id].m(); i++)
            {
              for(unsigned int j=0; j<coefficients[id].n(); j++)
              {
                out << std::setprecision(5) << std::scientific 
                  << coefficients[id](i,j) << " " ;
              } // j
              out << std::endl;
            } // i
          } // id
        } // if rank==0

      }


      /* unsigned int myRank = _comm.rank(); */
      /* unsigned int commSize = _comm.size(); */
      /* unsigned int coeffStart = std::min(_comm.rank(), _totalNCoeff-1); */

      /* if ( myRank == 0) */
      /*   out << "#" << std::string(75,'=') << std::endl; */

      /* for (unsigned int i=0; i < solutionNames.size(); i++) */
      /* { */
      /*   std::string id = solutionNames[i]; */

      /*   if (myRank == 0) */
      /*   { */
      /*     out << "#" << std::string(75,'-') << std::endl; */
      /*     out << "#" << "\t Solution: " << id << std::endl; */
      /*     out << "#" << std::string(75,'-') << std::endl; */
      /*   } */

      /*   for(unsigned int c=0; c<this->_totalNCoeff; c++) */
      /*   { */
      /*     unsigned int myC = (c-coeffStart)/commSize; */
      /*     std::vector<double> myCoeff(_solSize[id],0.); */
      /*     libMesh::Parallel::MessageTag tag(c); */
      /*     libMesh::Parallel::Status stat; */

      /*     if (c%commSize == myRank) */
      /*     { */
      /*       if (myRank == 0) */
      /*       { */
      /*         for(unsigned int comp=0; comp < (_coefficients[id])[myC].size(); comp++) */
      /*           out << std::setprecision(5) << std::scientific */ 
      /*             << (_coefficients[id])[myC](comp) << " " ; */
      /*         out << std::endl; */
      /*       } */
      /*       else */
      /*       { */
      /*         for(unsigned int comp=0; comp<myCoeff.size(); comp++) */
      /*           myCoeff[comp]= _coefficients[id][myC](comp) ; */ 
      /*         _comm.send(0,myCoeff,tag); */

      /*       } */
      /*     } */
      /*     else */
      /*     { */
      /*       if (myRank == 0) */
      /*       { */
      /*         stat=_comm.receive( c%commSize, myCoeff, tag); */
      /*         for(unsigned int comp=0; comp < myCoeff.size(); comp++) */
      /*           out << std::setprecision(5) << std::scientific */ 
      /*             << myCoeff[comp] << " " ; */
      /*         out << std::endl; */
      /*       } */
      /*     } */


      /*   } */
      /* } */


    }

/********************************************//**
 * \brief print coefficient values
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModelBase<T_S,T_P>::printCoefficients( 
        std::string solutionName,
        std::ostream& out ) 
    {
      std::vector<std::string> solutionsToPrint(1,solutionName);
      printCoefficients(solutionsToPrint, out);
    }

/********************************************//**
 * \brief print coefficient values
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModelBase<T_S,T_P>::printCoefficients( 
        std::ostream& out ) 
    {
      std::vector< std::string > solutionsToPrint ;
      
      std::set< std::string >::iterator id;
      for (id=_solutionNames.begin(); id!=_solutionNames.end(); id++)
        solutionsToPrint.push_back( *id ) ;
      printCoefficients(solutionsToPrint, out);
    }


/********************************************//**
 * \brief basic evaluation routine, based on more general (virtual) evaluate
 * routine
 ***********************************************/
  template<class T_S, class T_P> 
    T_P SurrogateModelBase<T_S,T_P>::evaluate( 
        std::string solutionName,  ///< solution to return
        T_S& parameterValues,    ///< parameter values to evaluate*/
        bool saveLocal  /**< save solution locally after evaluation*/
        ) const
    {
      std::set< std::string > solutionsToGet;
      solutionsToGet.insert(solutionName);
      std::map< std::string, T_P > solutionVectors
        = this->evaluate( solutionsToGet, parameterValues, saveLocal ) ;

      return solutionVectors[solutionName];
    }

/********************************************//**
 * \brief basic evaluation routine, based on more general (virtual) evaluate
 * routine
 ***********************************************/
  template<class T_S, class T_P> 
    std::map<std::string, T_P> SurrogateModelBase<T_S,T_P>::evaluate( 
        T_S& parameterValues,     ///< parameter values to evaluate*/
        bool saveLocal /**< save solution locally after evaluation*/
        ) const
    {
      std::set< std::string > solutionsToGet = _solutionNames;
      
      return evaluate( solutionsToGet, parameterValues, saveLocal ) ;    
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModelBase<T_S,T_P>::sample(  
          std::string solutionName, unsigned int N, std::vector<T_P>& sampleVec
        )
    {
      //clear out any old samples
      sampleVec.clear();

      // make sure solutionName is available
      agnos_assert( (_solutionNames.count( solutionName )) ) ;

      // variables needed for gsl randum number generator
      std::vector<const gsl_rng_type *> T;
      std::vector<gsl_rng *> r;
      gsl_rng_env_setup();

      // loop through parameters and set up rng
      for(unsigned int p=0; p<_parameters.size(); p++)
      {
        // initialize rng for this parameter type
        T.push_back( gsl_rng_default );
        r.push_back( gsl_rng_alloc(T.back()) );
      } // end 

      // sample generation loop
      for(unsigned int i=0; i<N; i++)
      {
        std::vector<double> s;
        for(unsigned int p=0; p<_parameters.size(); p++)
        {
          double scaledSample;
          // check type 
          // UNIFROM
          switch( ParameterType(_parameters[p]->type()) )
          {
            case UNIFORM:
              scaledSample = _parameters[p]->min() 
                + (_parameters[p]->max() - _parameters[p]->min() ) 
                * gsl_rng_uniform(r[p]) ;
              s.push_back( scaledSample  );
              break;
            case CONSTANT:
              s.push_back( _parameters[p]->min() ) ;
              break;
          }
        }

        // store as template vector type
        T_S paramValues(s);

        // evaluate model
        sampleVec.push_back( 
            this->evaluate( solutionName, paramValues, true) 
            ) ;


      } // end loop over samples

      // free up memory
      for(unsigned int p=0; p<_parameters.size(); p++)
        gsl_rng_free(r[p]);

    }


/********************************************//**
 * \brief calculate l2norm for single function
 * calls virtual l2Norm( solutionNames) routine
 ***********************************************/
  template<class T_S, class T_P> 
    T_P SurrogateModelBase<T_S,T_P>::l2Norm( 
        std::string solutionName  ///< solution to return
        ) 
    {
      std::set< std::string > solutionsToGet;
      solutionsToGet.insert(solutionName);

        std::map< std::string, T_P > solutionVectors
          = this->l2Norm( solutionsToGet ) ;

      return solutionVectors[solutionName];
    }

/********************************************//**
 * \brief calculate l2norm for all functions
 * calls virtual l2Norm( solutionNames) routine
 ***********************************************/
  template<class T_S, class T_P> 
    std::map<std::string, T_P> SurrogateModelBase<T_S,T_P>::l2Norm( ) 
    {
      std::set< std::string > solutionsToGet  = _solutionNames;
      return l2Norm( solutionsToGet) ;    
    }
  
  /********************************************//**
   * \brief get mean of surrogate model
   *      returns first coeff vector
   ***********************************************/
  template<class T_S,class T_P>
     std::map< std::string, T_P > SurrogateModelBase<T_S,T_P>::mean( )
    {
      std::map< std::string, T_P > meanCoefficients;

      // some initialization
      unsigned int myRank = _comm.rank();

      // mean is the first coefficient
      unsigned int coeff = 0 ;

      // some initialization
      libMesh::Parallel::MessageTag tag(coeff);
      libMesh::Parallel::Status stat;

      std::set<std::string>::iterator id = _solutionNames.begin();
      for (; id!=_solutionNames.end(); ++id)
      {
        unsigned int solSize = _solSize[*id];
        std::vector<double> solutionCoeff( solSize,0. );

        if (myRank == 0) // 0 coeff is on rank 0
        {
          // set value of mean coefficient
          for(unsigned int comp=0; comp<solSize; comp++)
            solutionCoeff[comp]= (*_coefficients[*id])(0,comp) ; 
        }

        // if I don't have coeff 0 receive it if I do send it
        _comm.broadcast(solutionCoeff);

        // convert to T_P type
        T_P coeffVector(solutionCoeff.size());
        for(unsigned int i=0; i<coeffVector.size(); i++)
          coeffVector(i) = solutionCoeff[i] ;

        meanCoefficients.insert(
            std::pair< std::string, T_P >( *id, coeffVector )
            );
      }

      return meanCoefficients;
    }
  

  /********************************************//**
   * \brief calculate l2 norm by integration
   ***********************************************/
  template<class T_S, class T_P>
    std::map< std::string, T_P> SurrogateModelBase<T_S,T_P>::l2Norm(
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

      if(this->_groupRank==0)
      {
        //---------
        //get integration points for higher order quad rule
        std::vector<unsigned int> integrationOrder = this->_order;
        for(unsigned int i=0; i<integrationOrder.size();i++)
          integrationOrder[i] += 0;

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
            = this->evaluate(solutionNames,paramValues,true);
          for (id=solutionNames.begin(); id != solutionNames.end(); id++)
          {

            for(unsigned int j=0; j<l2norm[*id].size();j++)
              l2norm[*id](j) += pointValue[*id](j) * pointValue[*id](j) 
                * quadWeights[i] ; 
          }

        } // end loop over quad points

      } // end groupRank == 0

      
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
    double SurrogateModelBase<T_S,T_P>::l2NormDifference(
        SurrogateModelBase<T_S,T_P>& comparisonModel,
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

      if(this->_groupRank==0)
      {
        //---------
        //get integration points for higher order quad rule
        std::vector<unsigned int> integrationOrder = this->_order;
        std::vector<unsigned int> comparisonOrder 
          = comparisonModel.getExpansionOrder() ;
        for(unsigned int i=0; i<integrationOrder.size();i++)
          integrationOrder[i] 
            = std::max(integrationOrder[i],comparisonOrder[i]) 
            + 0;

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
          /* std::cout << " test: base model eval " << std::endl; */
          diffVec -= comparisonModel.evaluate(solutionName,paramValues);
          /* std::cout << " test: comparison model eval " << std::endl; */


          l2norm += diffVec.dot( diffVec )  * quadWeights[i] ; 
          /* l2norm += 1  * quadWeights[i] ; */ 

        } // end loop over quad points

      } // end of if groupRank == 0

      
      // take square root
      l2norm =  std::sqrt( l2norm );


      return l2norm;

    }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
    double SurrogateModelBase<T_S,T_P>::l2NormDifference(
        std::string comparisonName,
        std::string solutionName )
    {

      // initialize to zero
      double l2norm = 0;

      // make sure both solutionName and comparisonName are present 
      // evaluate at integration points

      if (this->_solutionNames.count(solutionName) == 0) 
      {
        std::cerr << "\n\t ERROR: requested solution name not present in "
          << " base model of l2NormDifference( ... ). \n "
          << std::endl;
        exit(1);
      }
      if (this->_solutionNames.count(comparisonName) == 0 )
      {
        std::cerr << "\n\t ERROR: requested comparison name not present in "
          << " base model of l2NormDifference( ... ). \n "
          << std::endl;
        exit(1);
      }

      if(this->_groupRank==0)
      {
        //---------
        //get integration points for higher order quad rule
        std::vector<unsigned int> integrationOrder = this->_order;

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
          /* std::cout << " test: base model eval " << std::endl; */
          diffVec -= this->evaluate(comparisonName,paramValues);
          /* std::cout << " test: comparison model eval " << std::endl; */


          l2norm += diffVec.dot( diffVec )  * quadWeights[i] ; 
          /* l2norm += 1  * quadWeights[i] ; */ 

        } // end loop over quad points

      } // end of if groupRank == 0

      
      // take square root
      l2norm =  std::sqrt( l2norm );


      return l2norm;

    }

  /********************************************//**
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
    double SurrogateModelBase<T_S,T_P>::evaluateError(
        std::string solutionName )
    {

      // initialize to zero
      double error = 0;

      // make sure solutionName is present in surrogate and physics model
      if (this->_solutionNames.count(solutionName) == 0) 
      {
        std::cerr << "\n\t ERROR: requested solution name not present in "
          << " surrogate in evaluateError( ... ). \n "
          << std::endl;
        exit(1);
      }
      if (this->_physics->getAvailableSolutions().count(solutionName) == 0 )
      {
        std::cerr << "\n\t ERROR: requested solution name not present in "
          << " physics in evaluateError( ... ). \n "
          << std::endl;
        exit(1);
      }

      if(this->_groupRank==0)
      {
        //---------
        //get integration points for higher order quad rule
        std::vector<unsigned int> integrationOrder = this->_order;
        for (unsigned int i=0; i<integrationOrder.size(); i++)
          integrationOrder[i] *= 2.;

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

          // have to define a 1 element set with solution name in order to call
          // compute function from base PhysicsModel
          std::set<std::string> computeSolutions;
          computeSolutions.insert(solutionName);
          std::map<std::string, T_P> solutionVectors;

          // call physics compute to evaluate actual value of solution at this
          // parameter value
          this->_physics->compute( computeSolutions, paramValues, solutionVectors );

          T_P diffVec = solutionVectors[solutionName];
          diffVec -= this->evaluate(solutionName,paramValues);

          error += diffVec.dot( diffVec )  * quadWeights[i] ; 

        } // end loop over quad points

      } // end of if groupRank == 0

      
      // take square root
      error =  std::sqrt( error );


      return error;

    }

  /********************************************//**
   * \brief 
   ***********************************************/
    template<class T_S, class T_P> 
      void SurrogateModelBase<T_S,T_P>::printIndexSet( std::ostream& out ) const
      {
        if (this->_groupRank == 0)
        {
          out << std::endl;
      out << "#----------------------------------------------------" <<
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
    void SurrogateModelBase<T_S,T_P>::prettyPrintIndexSet( ) const
    {
      if (this->_groupRank == 0)
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



  template class 
    SurrogateModelBase<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;
} // namespace AGNOS
