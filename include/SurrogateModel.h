
// should be able to build up multi-element approaches based on basic surrogate
// definition I think. Even if we use libmesh?

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
          PhysicsModel<T_S,T_P>*            physics,
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
      SurrogateModel( 
          const SurrogateModel<T_S,T_P>* primarySurrogate, 
          unsigned int increaseOrder = 0,
          unsigned int multiplyOrder = 1,
          std::set<std::string> evaluateSolutions = std::vector<std::string>(),
          std::set<std::string> computeSolutions = std::vector<std::string>()
          );

      /** default destructor */
      virtual ~SurrogateModel( ); 

      /** build the surrogate model construction */
      virtual void build( ) = 0; 

      /** evaluate surrogate model at give parameterValues and return
       * solutionNames */
      virtual std::map<std::string, T_P> evaluate( 
          std::set<std::string> solutionNames,  ///< solution to return
          T_S& parameterValues /**< parameter values to evaluate*/
          ) const = 0;
      /** evaluate surrogate model at give parameterValues and return
       * solutionName */
      T_P evaluate( 
          std::string solutionName,  ///< solution to return
          T_S& parameterValues     ///< parameter values to evaluate*/
          ) const ;
      /** evaluate surrogate model at give parameterValues and return
       * all solutions*/
      std::map<std::string, T_P> evaluate( 
          T_S& parameterValues     ///< parameter values to evaluate*/
          ) const ;

      /** Refine the surrogate model. Must be definied in derived classes. */
      virtual void refine( ) = 0;

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
      void setParameters( std::vector<Parameter*> parameters );
      /** return reference to parameters object */
      std::vector<Parameter*>   getParameters( ) const
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
      const std::map< std::string, std::vector<T_P> >   getLocalCoefficients() const
      { return _coefficients; }
      /** reference to all coefficients */
      const std::map< std::string, std::vector<T_P> >   getCoefficients() ;

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

      /** expansion order used to construct surrogateModel */
      const std::vector<unsigned int> getExpansionOrder( ) const
      { return _order; }

      /** reference to communicator  */
      const Communicator& getComm( ) const 
      { return _comm; }
      /** reference to physics pointer */
      PhysicsModel<T_S,T_P>* getPhysics( ) const
      { return _physics; }

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
      PhysicsModel<T_S,T_P>* _physics;
      
      /** expansion order */
      std::vector<unsigned int>                           _order;  
      /** coefficients vectors (local)*/
      std::map< std::string, std::vector<T_P> >           _coefficients;
      /** total number of coefficients on all vectors */
      unsigned int                                        _totalNCoeff;

      /** dimension of coefficient vectors for each solution name */
      std::map< std::string, unsigned int>                _solSize;

      /** reference to Parameter object */
      std::vector<Parameter*>                             _parameters;
      /** Parameter dimension */
      unsigned int                                        _dimension;

      /** Set of solution names that the surrogate model computes */
      std::set<std::string>                               _solutionNames;

      /** Set of solution names that the surrogate model usese primary surrogate
       * for. Only meaningful in secondary constructor case */
      std::set<std::string>                               _evalNames;

      /** Primary surrogate to use in evaluation for secondary surrogate
       * construciton */
      const SurrogateModel<T_S,T_P>*                      _evalSurrogate ;

      /** Data structure to hold evalSurrogate evaluations, to be used in
       * surrogate construction */
      std::vector<std::map< std::string,T_P> > _primalEvaluations;
      

  }; //SurrogateModel class

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
        const Communicator&               comm,
        PhysicsModel<T_S,T_P>*            physics,
        const std::vector<Parameter*>     parameters,
        const std::vector<unsigned int>&          order
        )
      : 
        _comm(comm),
        _physics(physics), 
        _parameters(parameters),
        _dimension( parameters.size() ), 
        _order(order)
    {
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
    }

/********************************************//*
 * \brief Secondary constructor
 ***********************************************/
  template<class T_S, class T_P>
    SurrogateModel<T_S,T_P>::SurrogateModel( 
          const SurrogateModel<T_S,T_P>* primarySurrogate, 
          unsigned int increaseOrder ,
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
      // augment order appropriately
      _order = primarySurrogate->getExpansionOrder() ;
      for (unsigned int i=0; i<_order.size();i++)
      {
        _order[i] += increaseOrder ;
        _order[i] *= multiplyOrder ;
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
    }

/********************************************//**
 * \brief set parameters
 ***********************************************/
  template<class T_S, class T_P>
    void SurrogateModel<T_S,T_P>::setParameters( 
        std::vector<Parameter*> parameters)
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
    const std::map< std::string, std::vector<T_P> >
    SurrogateModel<T_S,T_P>::getCoefficients( ) 
    {
     
      std::map< std::string, std::vector<T_P> > allCoefficients;

      unsigned int myRank = _comm.rank();
      unsigned int commSize = _comm.size();
      unsigned int coeffStart = std::min(_comm.rank(), _totalNCoeff-1);

      // loop over all solution names
      std::set<std::string>::iterator id = _solutionNames.begin();
      for (; id!=_solutionNames.end(); ++id)
      {
        unsigned int solSize = _solSize[*id];
        std::vector<T_P> solutionCoeff( _totalNCoeff, T_P(solSize) );

        for(unsigned int c=0; c<this->_totalNCoeff; c++)
        {
          unsigned int myC = (c-coeffStart)/commSize;
          std::vector<double> myCoeff(solSize,0.) ;
          libMesh::Parallel::MessageTag tag(c);
          libMesh::Parallel::Status stat;

          if (c%commSize == myRank)
          {
            if (myRank == 0)
            {
              // get and insert in total set
              solutionCoeff[c] = _coefficients[*id][myC] ;
            }
            else
            {
              // send my coeff to rank 0
              for(unsigned int comp=0; comp<myCoeff.size(); comp++)
                myCoeff[comp]= _coefficients[*id][myC](comp) ; 
              _comm.send(0,myCoeff,tag);
            }
          }
          else
          {
            if (myRank == 0)
            {
              // if im rank 0 receive all coefficients
              stat=_comm.receive( c%commSize, myCoeff, tag);

              // insert into global set
              for(unsigned int i=0; i<solSize; i++)
                solutionCoeff[c](i) = myCoeff[i];

            }
          }


        }


        allCoefficients.insert( 
            std::pair<std::string, std::vector<T_P> >(
              *id, solutionCoeff )
            );
      }

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
      unsigned int myRank = _comm.rank();
      unsigned int commSize = _comm.size();
      unsigned int coeffStart = std::min(_comm.rank(), _totalNCoeff-1);

      if ( myRank == 0)
        out << "#" << std::string(75,'=') << std::endl;

      for (unsigned int i=0; i < solutionNames.size(); i++)
      {
        std::string id = solutionNames[i];

        if (myRank == 0)
        {
          out << "#" << std::string(75,'-') << std::endl;
          out << "#" << "\t Solution: " << id << std::endl;
          out << "#" << std::string(75,'-') << std::endl;
        }

        for(unsigned int c=0; c<this->_totalNCoeff; c++)
        {
          unsigned int myC = (c-coeffStart)/commSize;
          std::vector<double> myCoeff(_solSize[id],0.);
          libMesh::Parallel::MessageTag tag(c);
          libMesh::Parallel::Status stat;

          if (c%commSize == myRank)
          {
            if (myRank == 0)
            {
              for(unsigned int comp=0; comp < (_coefficients[id])[myC].size(); comp++)
                out << std::setprecision(5) << std::scientific 
                  << (_coefficients[id])[myC](comp) << " " ;
              out << std::endl;
            }
            else
            {
              for(unsigned int comp=0; comp<myCoeff.size(); comp++)
                myCoeff[comp]= _coefficients[id][myC](comp) ; 
              _comm.send(0,myCoeff,tag);

            }
          }
          else
          {
            if (myRank == 0)
            {
              stat=_comm.receive( c%commSize, myCoeff, tag);
              for(unsigned int comp=0; comp < myCoeff.size(); comp++)
                out << std::setprecision(5) << std::scientific 
                  << myCoeff[comp] << " " ;
              out << std::endl;
            }
          }


        }
      }


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
        T_S& parameterValues     ///< parameter values to evaluate*/
        ) const
    {
      std::set< std::string > solutionsToGet;
      solutionsToGet.insert(solutionName);
      std::map< std::string, T_P > solutionVectors
        = this->evaluate( solutionsToGet, parameterValues ) ;

      return solutionVectors[solutionName];
    }

/********************************************//**
 * \brief basic evaluation routine, based on more general (virtual) evaluate
 * routine
 ***********************************************/
  template<class T_S, class T_P> 
    std::map<std::string, T_P> SurrogateModel<T_S,T_P>::evaluate( 
        T_S& parameterValues     ///< parameter values to evaluate*/
        ) const
    {
      std::set< std::string > solutionsToGet = _solutionNames;
      
      return evaluate( solutionsToGet, parameterValues ) ;    
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
      unsigned int commSize = _comm.size();
      unsigned int coeffStart = std::min(_comm.rank(), _totalNCoeff-1);

      // mean is the first coefficient
      unsigned int coeff = 0 ;

      // some initialization
      unsigned int myC = (coeff-coeffStart)/commSize;
      libMesh::Parallel::MessageTag tag(coeff);
      libMesh::Parallel::Status stat;

      std::set<std::string>::iterator id = _solutionNames.begin();
      for (; id!=_solutionNames.end(); ++id)
      {
        std::vector<double> myCoeff(_solSize[*id],0.) ;
        T_P solutionCoeff( _solSize[*id] );

        if (coeff%commSize == myRank)
        {
          // set value of mean coefficient
          for(unsigned int comp=0; comp<myCoeff.size(); comp++)
            myCoeff[comp]= _coefficients[*id][myC](comp) ; 
        }

        // if I don't have coeff 0 receive it if I do send it
        _comm.broadcast(myCoeff, coeff%commSize);

        // insert into global set
        for(unsigned int i=0; i<solutionCoeff.size(); i++)
          solutionCoeff(i) = myCoeff[i];

        meanCoefficients.insert(
            std::pair< std::string, T_P>( *id, solutionCoeff )
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
