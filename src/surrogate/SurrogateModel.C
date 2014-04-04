
#include "SurrogateModel.h"

namespace AGNOS
{


/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        const Communicator&               comm,
        std::shared_ptr<PhysicsModel<T_S,T_P> >           physics,
        const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
        const std::vector<unsigned int>&          order,
        std::set<std::string> computeSolutions 
        )
      : 
        _comm(comm),
        _physics(physics), 
        _parameters(parameters),
        _dimension( parameters.size() ),
        _order(order),
        _solutionNames(computeSolutions)
    {
      int globalRank,globalSize;
      MPI_Comm_rank(MPI_COMM_WORLD,&globalRank);
      MPI_Comm_size(MPI_COMM_WORLD,&globalSize);
      _physicsGroup = globalRank / ( _physics->comm().size() ) ;
      _nPhysicsGroups = globalSize / ( _physics->comm().size() ) ;
      MPI_Comm_rank(_physics->comm().get(),&_groupRank);

      assert( _groupRank == (globalRank % _physics->comm().size() ) );


      //check against available physics solutions
      std::set<std::string> physicsSolutions = physics->getAvailableSolutions();
      if (_solutionNames.empty())
        _solutionNames = physicsSolutions ;
      else
      {
        std::set<std::string>::iterator solName = _solutionNames.begin() ;
        for( ; solName != _solutionNames.end(); solName++)
          if ( !physicsSolutions.count(*solName) ) 
          {
            std::cerr << std::endl;
            std::cerr 
              << "ERROR: requested solution "
              << *solName << " " 
              << " is not available in this PhysicsModel class. \n" ;
            std::cerr << std::endl;
            std::abort();
          }
      }
      

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

/********************************************//*
 * \brief Secondary constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
          std::shared_ptr<SurrogateModel<T_S,T_P> > primarySurrogate, 
          std::vector<unsigned int> increaseOrder ,
          unsigned int multiplyOrder ,
          std::set<std::string> evaluateSolutions , 
          std::set<std::string> computeSolutions 
        )
      : 
        _comm( primarySurrogate->getComm() ),
        _physics( primarySurrogate->getPhysics() ), 
        _parameters( primarySurrogate->getParameters() ),
        _dimension( _parameters.size() ),
        _evalSurrogate( primarySurrogate )
    {
      int globalRank,globalSize;
      MPI_Comm_rank(MPI_COMM_WORLD,&globalRank);
      MPI_Comm_size(MPI_COMM_WORLD,&globalSize);
      _physicsGroup = globalRank / ( _physics->comm().size() ) ;
      _nPhysicsGroups = globalSize / ( _physics->comm().size() ) ;
      MPI_Comm_rank(_physics->comm().get(),&_groupRank);

      assert( _groupRank == (globalRank % _physics->comm().size() ) );


      // augment order appropriately
      _order = primarySurrogate->getExpansionOrder() ;
      if (increaseOrder.empty())
        _increaseOrder = std::vector<unsigned int>(_dimension,0) ;
      else
        _increaseOrder = increaseOrder;

      agnos_assert(_increaseOrder.size() == _dimension);

      // check/reset increase order for CONSTANT parameters
        for( unsigned int i=0; i<_parameters.size(); i++ )
          if ( _parameters[i]->type() == CONSTANT )
            if (_increaseOrder[i] != 0)
            {
              std::cout 
                << "WARNING: forcing CONSTANT parameter (" 
                << i 
                << ") increaseOrder to 0 \n" ;
              _increaseOrder[i] = 0;
            }

      _multiplyOrder = multiplyOrder ;
      for (unsigned int i=0; i<_order.size();i++)
      {
        _order[i] += _increaseOrder[i] ;
        _order[i] *= _multiplyOrder ;
      }

      // if no output solutions are provided assume we want all 
      if (computeSolutions.size() == 0)
        _solutionNames = primarySurrogate->getSolutionNames();
      else
        _solutionNames  = computeSolutions ;

      // if no evaluateSolutions are provided assume we use all
      if (evaluateSolutions.size() == 0)
        _evalNames = _solutionNames ;
      else
        _evalNames = evaluateSolutions ;
      

    }


/********************************************//**
 * \brief Default destructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::~SurrogateModel( )
    {
      _coefficients.clear();
    }

/********************************************//**
 * \brief set parameters
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModel<T_S,T_P>::setParameters( 
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
    SurrogateModel<T_S,T_P>::getCoefficients( ) 
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

          unsigned int solSize = _solSize[*id];
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
                solCoefficients(i,j) = (*this->_coefficients[*id])(i,j)  ;
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
                  cbuf.push_back( (*this->_coefficients[*id])(i,j) )  ;
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
    void SurrogateModel<T_S,T_P>::printCoefficients( 
        std::vector<std::string> solutionNames,
        std::ostream& out ) 
    {

      if ( this->_groupRank == 0)
      {
        std::map<std::string,LocalMatrix > coefficients =
          this->getCoefficients();

        if ( this->_physicsGroup == 0)
        {
          out << "#" << std::string(75,'=') << std::endl;

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
    void SurrogateModel<T_S,T_P>::printCoefficients( 
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
    void SurrogateModel<T_S,T_P>::printCoefficients( 
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
    T_P SurrogateModel<T_S,T_P>::evaluate( 
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
    std::map<std::string, T_P> SurrogateModel<T_S,T_P>::evaluate( 
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
    void SurrogateModel<T_S,T_P>::refine( )
    {
      std::vector<unsigned int> increase(this->_dimension, 1);
      refine( increase );
    }


/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModel<T_S,T_P>::refine( 
        const std::vector<unsigned int>& increase 
        )
    {

      // verify correct dimension of increase is given
      agnos_assert( (increase.size() == _parameters.size()) );

      // if this is a primary surrogate incease order
      if ( _evalSurrogate == NULL )
      {
        for( unsigned int i=0; i<_parameters.size(); i++ )
        {
          // increase order
          this->_order[i] += increase[i] ;

          // if its a constant parameter and we raised the order reset it back
          // to zero
          if ( _parameters[i]->type() == CONSTANT )
            if (this->_order[i] != 0)
            {
              std::cout 
                << "WARNING: forcing CONSTANT parameter (" 
                << i 
                << ") order to 0th order \n" ;
              this->_order[i] = 0;
            } //end if increase != 0
          //end if paramType = CONSTANT
          
        } //end for each parameter
      }
      // otherwise refine based on primary surrogate
      else
      {
        _order = _evalSurrogate->getExpansionOrder();
        for (unsigned int i=0; i<_order.size();i++)
        {
          _order[i] += _increaseOrder[i] ;
          _order[i] *= _multiplyOrder ;
        }

      }

      /* this->initialize(); */
      /* this->build(); */
    }

/********************************************//**
 * \brief calculate l2norm for single function
 * calls virtual l2Norm( solutionNames) routine
 ***********************************************/
  template<class T_S, class T_P> 
    T_P SurrogateModel<T_S,T_P>::l2Norm( 
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
    std::map<std::string, T_P> SurrogateModel<T_S,T_P>::l2Norm( ) 
    {
      std::set< std::string > solutionsToGet  = _solutionNames;
      return l2Norm( solutionsToGet) ;    
    }
  
  /********************************************//**
   * \brief get mean of surrogate model
   *      returns first coeff vector
   ***********************************************/
  template<class T_S,class T_P>
     std::map< std::string, T_P > SurrogateModel<T_S,T_P>::mean( )
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
            solutionCoeff[comp]= (*_coefficients[*id])(comp,0) ; 
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
   * \brief 
   ***********************************************/
  template<class T_S,class T_P>
    double SurrogateModel<T_S,T_P>::l2NormDifference(
        SurrogateModel<T_S,T_P>& comparisonModel,
        std::string solutionName )
    {
      std::cerr << "\n\t ERROR: l2NormDifference( ... )"
        << " is not implemented in derived class\n" 
        << std::endl;
      exit(1);
    }

  template class 
    SurrogateModel<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;


}
