
#ifndef SURROGATE_MODEL_H 
#define SURROGATE_MODEL_H 

#include "agnosDefines.h"
#include "Parameter.h"
#include "PhysicsModel.h"

namespace AGNOS
{

  enum SurrogateModelType{
    PSEUDO_SPECTRAL_TENSOR_PRODUCT=0,
    PSEUDO_SPECTRAL_SPARSE_GRID,
    PSEUDO_SPECTRAL_MONTE_CARLO,
    COLLOCATION };

  /********************************************//**
   * \brief Base surrogate model class
   *
   * Allows for derivation of Collocation and Pseudospectral surrogate models.
   * Function to construct SurrogateModel for must me defined by providing a
   * PhysicsFunction object. 
   ***********************************************/
  template<class T_S, class T_P>
  class SurrogateModel
  {


    public: 

      /** Constructor */
      SurrogateModel(
          const Communicator&               comm,
          std::shared_ptr<PhysicsModel<T_S,T_P> >            physics,
          const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
          const std::vector<unsigned int>&  order
          );

      /** Secondary Constructor. 
       *  ***DO NOT MISTAKE FOR A COPY CONSTRUCTOR ***
       *  Intended use is for constructing a secondary surrogate model using the
       *  primary model as an evaluating object in the build routine. 
       *  If additional inputs are defined it will
       * construct a new surrogate increasing the order and using
       * primarySurrogate to perform evaluations in the constructions */
      SurrogateModel( 
          std::shared_ptr<SurrogateModel<T_S,T_P> > primarySurrogate, 
          std::vector<unsigned int> increaseOrder = std::vector<std::string>(),
          unsigned int multiplyOrder = 1,
          std::set<std::string> evaluateSolutions = std::vector<std::string>(),
          std::set<std::string> computeSolutions = std::vector<std::string>()
          );

      /** default destructor */
      virtual ~SurrogateModel( ); 

      /** Initialization routine */
      virtual void initialize( ) = 0 ;

      /** build the surrogate model construction */
      virtual void build( ) = 0; 

      /** evaluate surrogate model at give parameterValues and return
       * solutionNames */
      virtual std::map<std::string, T_P> evaluate( 
          std::set<std::string> solutionNames,  ///< solution to return
          T_S& parameterValues, /**< parameter values to evaluate*/
          bool saveLocal = true /**< save solution locally after evaluation*/
          ) const = 0;
      /** evaluate surrogate model at give parameterValues and return
       * solutionName */
      T_P evaluate( 
          std::string solutionName,  ///< solution to return
          T_S& parameterValues,    ///< parameter values to evaluate*/
          bool saveLocal = true /**< save solution locally after evaluation*/
          ) const ;
      /** evaluate surrogate model at give parameterValues and return
       * all solutions*/
      std::map<std::string, T_P> evaluate( 
          T_S& parameterValues,    ///< parameter values to evaluate*/
          bool saveLocal = true /**< save solution locally after evaluation*/
          ) const ;

      /** Refine the surrogate model. Must be definied in derived classes. */
      virtual void refine( 
          const std::vector<unsigned int>& increase 
            = std::vector<unsigned int>() 
          ) ;

      /** calculate mean */
      std::map< std::string, T_P > mean( ) ;

      /** calculate the L2 norm over parameter space */
      virtual std::map<std::string, T_P> l2Norm( 
          std::set<std::string> solutionNames  ///< solution to return
          ) = 0;
      /** calculate the L2 norm over parameter space */
      T_P l2Norm( 
          std::string solutionName  ///< solution to return
          ) ;
      /** calculate the L2 norm over parameter space */
      std::map<std::string, T_P> l2Norm( ) ;

      /** Compute the norm of the difference between this and comaprisonModel,
       * for given solutionName */
      virtual double l2NormDifference( 
          SurrogateModel<T_S,T_P>& comparisonModel,
          std::string solutionName
          ) ;


      /** set parameters object  */
      void setParameters( std::vector<std::shared_ptr<AGNOS::Parameter> >& parameters );
      /** return reference to parameters object */
      const std::vector<std::shared_ptr<AGNOS::Parameter> >   getParameters( ) const
      { return _parameters; }

      /** print the coefficient vectors */
      void printCoefficients( 
        std::vector<std::string> solutionNames,
        std::ostream& out ) ;
      /** print the coefficient vectors */
      void printCoefficients( std::string solutionName, std::ostream& out ) ;
      /** print the coefficient vectors */
      void printCoefficients( std::ostream& out ) ;

      /** reference to locally stored coefficients */
      /* const std::map< std::string, LocalMatrix>   getLocalCoefficients() const */
      /* { return _coefficients; } */
      /** reference to all coefficients */
      const std::map< std::string, LocalMatrix>   getCoefficients() ;

      /** print integration weights */
      virtual void printIntegrationWeights( std::ostream& out ) const = 0;
      /** print integration points */
      virtual void printIntegrationPoints( std::ostream& out ) const = 0;
      /** print index set */
      virtual void printIndexSet( std::ostream& out ) const = 0;
      /** print integration weights in table format*/
      virtual void prettyPrintIntegrationWeights( ) const = 0;
      /** print integration weights in table format*/
      virtual void prettyPrintIntegrationPoints( ) const = 0;
      /** print integration weights in table format*/
      virtual void prettyPrintIndexSet( ) const = 0;

      /** Return index set for this surrogate */
      virtual const std::vector< std::vector<unsigned int> > indexSet() const = 0 ;

      /** expansion order used to construct surrogateModel */
      const std::vector<unsigned int> getExpansionOrder( ) const
      { return _order; }

      /** reference to communicator  */
      const Communicator& getComm( ) const 
      { return _comm; }
      /** reference to physics pointer */
      std::shared_ptr<PhysicsModel<T_S,T_P> > getPhysics( ) const
      { return _physics; }
      /** set reference to physics pointer */
      void setPhysics( std::shared_ptr<PhysicsModel<T_S,T_P> > physics ) 
      { _physics = physics; }
      /** reference to physicsGroup*/
      const int physicsGroup() const
      { return _physicsGroup; }
      /** reference to number of physics groups*/
      const int nPhysicsGroups() const
      { return _nPhysicsGroups; }
      /** reference to groupRank */
      const int groupRank() const
      { return _groupRank; }

      /** Parameter dimension */
      unsigned int dimension(){ return _dimension; }

      /** solution names this surrogateModel is built for */
      std::set<std::string> getSolutionNames( ) const
      { return _solutionNames; }

      /** reference to sol vector sizes */
      std::map<std::string, unsigned int> getSolSize( ) const
      { return _solSize; }

    protected: 
      /** reference to communicator */
      const Communicator& _comm;
      /** reference to underlying physics */
      std::shared_ptr<PhysicsModel<T_S,T_P> > _physics;
      /** reference to physics group number */
      int _physicsGroup;
      /** reference to number of physics groups */
      int _nPhysicsGroups;
      /** reference to groupRank (needed to deterimine if this is a master or
       * slave node) */
      int _groupRank;
      
      /** expansion order */
      std::vector<unsigned int> _order;  
      std::vector<unsigned int> _increaseOrder ;
      unsigned int              _multiplyOrder ;

      /** coefficients vectors */
      std::map< std::string, std::shared_ptr<DistMatrix> >                 _coefficients;
      /** total number of coefficients on all vectors */
      unsigned int                                        _totalNCoeff;
      /** indicies of coefficients corresponding to local process */
      std::vector<unsigned int> _coeffIndices;

      /** dimension of coefficient vectors for each solution name */
      std::map< std::string, unsigned int>                _solSize;

      /** reference to Parameter object */
      std::vector<std::shared_ptr<AGNOS::Parameter> >                   _parameters;
      /** Parameter dimension */
      unsigned int                                        _dimension;

      /** Set of solution names that the surrogate model computes */
      std::set<std::string>                               _solutionNames;

      /** Set of solution names that the surrogate model usese primary surrogate
       * for. Only meaningful in secondary constructor case */
      std::set<std::string>                               _evalNames;

      /** Primary surrogate to use in evaluation for secondary surrogate
       * construciton */
      std::shared_ptr<SurrogateModel<T_S,T_P> >      _evalSurrogate ;

      /** Data structure to hold evalSurrogate evaluations, to be used in
       * surrogate construction */
      std::vector< std::map< std::string,T_P> > _primaryEvaluations;
      

  }; //SurrogateModel class

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        const Communicator&               comm,
        std::shared_ptr<PhysicsModel<T_S,T_P> >           physics,
        const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
        const std::vector<unsigned int>&          order
        )
      : 
        _comm(comm),
        _physics(physics), 
        _parameters(parameters),
        _dimension( parameters.size() ),
        _order(order)
    {
      int globalRank,globalSize;
      MPI_Comm_rank(MPI_COMM_WORLD,&globalRank);
      MPI_Comm_size(MPI_COMM_WORLD,&globalSize);
      _physicsGroup = globalRank / ( _physics->comm().size() ) ;
      _nPhysicsGroups = globalSize / ( _physics->comm().size() ) ;
      MPI_Comm_rank(_physics->comm().get(),&_groupRank);

      assert( _groupRank == (globalRank % _physics->comm().size() ) );


      _solutionNames.clear();
      _solutionNames = physics->getSolutionNames();


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
    void SurrogateModel<T_S,T_P>::refine( 
        const std::vector<unsigned int>& increase 
        )
    {
      // if this is a primary surrogate incease order
      if ( _evalSurrogate == NULL )
      {
        if ( increase.empty() )
          for(unsigned int i=0; i<this->_dimension; i++)
            this->_order[i]++;
        else
          for(unsigned int i=0; i<this->_dimension; i++)
            this->_order[i] += increase[i] ;
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


}
#endif //SURROGATE_MODEL_H
