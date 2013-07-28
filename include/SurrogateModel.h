
// should be able to build up multi-element approaches based on basic surrogate
// definition I think. Even if we use libmesh?

#ifndef SURROGATE_MODEL_H 
#define SURROGATE_MODEL_H 

#include "agnosDefines.h"
#include "Parameter.h"

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
   * PhysicsModel object. 
   ***********************************************/
  template<class T_S>
  class SurrogateModel
  {


    public: 

      // single physics function constructors
      SurrogateModel(
          const Communicator&               comm,
          PhysicsModel<T_S>*                physics,
          const std::vector<Parameter*>     parameters,
          const unsigned int                order 
          );
      SurrogateModel(
          const Communicator&               comm,
          PhysicsModel<T_S>*                physics,
          const std::vector<Parameter*>     parameters,
          const std::vector<unsigned int>&  order
          );

      SurrogateModel( );           /**< Default constructor */
      virtual ~SurrogateModel( );  /**< Default destructor */

      // surrogate construction and evaluation
      virtual void build( ) = 0; 

      virtual std::map<std::string, std::vector<Number> > evaluate( 
          std::vector<std::string> solutionNames,  ///< solution to return
          T_S& parameterValues /**< parameter values to evaluate*/
          ) = 0;
      std::vector<Number> evaluate( 
          std::string solutionName,  ///< solution to return
          T_S& parameterValues     ///< parameter values to evaluate*/
          ) ;
      std::map<std::string, std::vector<Number> > evaluate( 
          T_S& parameterValues     ///< parameter values to evaluate*/
          ) ;

      virtual void refine( ) = 0;

      // calculate norms and statistics
      /* std::map< std::string, NumericVector<Number>  > mean( ) ; */
      /* // integration method depends on surr type */
      /* virtual std::map<std::string, NumericVector<Number> > l2Norm( */ 
      /*     std::vector<std::string> solutionNames  ///< solution to return */
      /*     ) = 0; */
      /* NumericVector<Number>  l2Norm( */ 
      /*     std::string solutionName  ///< solution to return */
      /*     ) ; */
      /* std::map<std::string, NumericVector<Number> > l2Norm( ) ; */

      /* // difference btw two surrogate models */
      /* virtual double l2NormDifference( */ 
      /*     SurrogateModel<T_S>& comparisonModel, */
      /*     std::string solutionName */
      /*     ) ; */


      // Manipulators
      void setParameters( std::vector<Parameter*> parameters );
      std::vector<Parameter*>   getParameters( ) const;

      /* void printCoefficients( */ 
      /*   std::vector<std::string> solutionNames, */
      /*   std::ostream& out ) ; */
      /* void printCoefficients( std::string solutionName, std::ostream& out ) ; */
      /* void printCoefficients( std::ostream& out ) ; */
      const std::map< std::string, std::vector< std::vector<Number> > >   getLocalCoefficients() const;
      /* const std::map< std::string, std::vector<NumericVector<Number> > >   getCoefficients() ; */

      virtual void printIntegrationWeights( std::ostream& out ) const = 0;
      virtual void printIntegrationPoints( std::ostream& out ) const = 0;
      virtual void printIndexSet( std::ostream& out ) const = 0;
      virtual void prettyPrintIntegrationWeights( ) const = 0;
      virtual void prettyPrintIntegrationPoints( ) const = 0;
      virtual void prettyPrintIndexSet( ) const = 0;

      std::vector<unsigned int> getExpansionOrder( ) const;

      std::set<std::string> getSolutionNames( ) const
        { return _solutionNames; };

    protected: 
      const Communicator& _comm;
      PhysicsModel<T_S>* _physics;
      
      std::vector<unsigned int>                           _order;  
      std::map< std::string, std::vector< std::vector<Number> > > _coefficients;
      unsigned int                                        _totalNCoeff;

      std::map< std::string, unsigned int>                _solSize;

      std::vector<Parameter*>                             _parameters;
      unsigned int                                        _dimension;

      std::set<std::string>                               _solutionNames;

      

  }; //SurrogateModel class

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S>
    SurrogateModel<T_S>::SurrogateModel( 
        const Communicator&               comm,
          PhysicsModel<T_S>*                physics,
        std::vector<Parameter*>                             parameters,
        unsigned int                                        order
        )
      : 
        _comm(comm),
        _physics(physics), 
        _parameters(parameters),
        _dimension( parameters.size() )
    {
      _order = std::vector<unsigned int>(_dimension,order);

      _solutionNames.clear();
      _solutionNames = physics->getSolutionNames();


    }

/********************************************//*
 * \brief Constructor
 ***********************************************/
  template<class T_S>
    SurrogateModel<T_S>::SurrogateModel( 
        const Communicator&               comm,
        PhysicsModel<T_S>*                physics,
        std::vector<Parameter*>                   parameters,
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
 * \brief Default constructor
 ***********************************************/
  template<class T_S>
    SurrogateModel<T_S>::SurrogateModel( )
    {
    }

/********************************************//**
 * \brief Default destructor
 ***********************************************/
  template<class T_S>
    SurrogateModel<T_S>::~SurrogateModel( )
    {
    }

/********************************************//**
 * \brief set parameters
 ***********************************************/
  template<class T_S>
    void SurrogateModel<T_S>::setParameters( 
        std::vector<Parameter*> parameters)
    {
      _parameters = parameters ;
      _dimension = parameters.size();
      return;
    }

/********************************************//**
 * \brief get number of parameters
 ***********************************************/
  template<class T_S>
    std::vector<Parameter*> SurrogateModel<T_S>::getParameters( ) const
    {
      return _parameters ;
    }

/********************************************//**
 * \brief Get the current expansion order
 ***********************************************/
  template<class T_S>
    std::vector<unsigned int> SurrogateModel<T_S>::getExpansionOrder( )
    const
  {
    return _order;
  }

/********************************************//**
 * \brief return coefficient values
 ***********************************************/
  template<class T_S> 
    const std::map< std::string, std::vector< std::vector<Number> > >
    SurrogateModel<T_S>::getLocalCoefficients( ) const
    {
      return _coefficients;
    }

/********************************************//**
 * \brief return all coefficients to rank 0, to save on communication we do not
 * broadcast to all processes
 ***********************************************/
  /* template<class T_S> */ 
  /*   const std::map< std::string, std::vector<NumericVector<Number>* > > */
  /*   SurrogateModel<T_S>::getCoefficients( ) */ 
  /*   { */
     
  /*     std::map< std::string, std::vector<NumericVector<Number>* > > allCoefficients; */

  /*     unsigned int myRank = _comm.rank(); */
  /*     unsigned int commSize = _comm.size(); */
  /*     unsigned int coeffStart = std::min(_comm.rank(), _totalNCoeff-1); */

  /*     // loop over all functions */
  /*     std::set<std::string>::iterator id = _solutionNames.begin(); */
  /*     for (; id!=_solutionNames.end(); ++id) */
  /*     { */
  /*       /1* unsigned int solSize = _solSize[*id]; *1/ */
  /*       std::vector<NumericVector<Number>* > solutionCoeff; */
  /*       solutionCoeff.resize(_totalNCoeff); */

  /*       for(unsigned int c=0; c<this->_totalNCoeff; c++) */
  /*       { */
  /*         unsigned int myC = (c-coeffStart)/commSize; */
  /*         std::vector<double> myCoeff(solSize,0.) ; */
  /*         libMesh::Parallel::MessageTag tag(c); */
  /*         libMesh::Parallel::Status stat; */

  /*         if (c%commSize == myRank) */
  /*         { */
  /*           if (myRank == 0) */
  /*           { */
  /*             // get and insert in total set */
  /*             solutionCoeff[c] = _coefficients[*id][myC] ; */
  /*           } */
  /*           else */
  /*           { */
  /*             // send my coeff to rank 0 */
  /*             for(unsigned int comp=0; comp<myCoeff.size(); comp++) */
  /*               myCoeff[comp]= _coefficients[*id][myC](comp) ; */ 
  /*             _comm.send(0,myCoeff,tag); */
  /*           } */
  /*         } */
  /*         else */
  /*         { */
  /*           if (myRank == 0) */
  /*           { */
  /*             // if im rank 0 receive all coefficients */
  /*             stat=_comm.receive( c%commSize, myCoeff, tag); */

  /*             // insert into global set */
  /*             for(unsigned int i=0; i<solSize; i++) */
  /*               solutionCoeff[c](i) = myCoeff[i]; */

  /*           } */
  /*         } */


  /*       } */


  /*       allCoefficients.insert( */ 
  /*           std::pair<std::string, std::vector<NumericVector<Number> > >( */
  /*             *id, solutionCoeff ) */
  /*           ); */
  /*     } */

  /*     return allCoefficients; */
  /*   } */

/********************************************//**
 * \brief print coefficient values
 ***********************************************/
  /* template<class T_S> */
  /*   void SurrogateModel<T_S>::printCoefficients( */ 
  /*       std::vector<std::string> solutionNames, */
  /*       std::ostream& out ) */ 
  /*   { */
  /*     unsigned int myRank = _comm.rank(); */
  /*     unsigned int commSize = _comm.size(); */
  /*     unsigned int coeffStart = std::min(_comm.rank(), _totalNCoeff-1); */

  /*     if ( myRank == 0) */
  /*       out << "#" << std::string(75,'=') << std::endl; */

  /*     for (unsigned int i=0; i < solutionNames.size(); i++) */
  /*     { */
  /*       std::string id = solutionNames[i]; */

  /*       if (myRank == 0) */
  /*       { */
  /*         out << "#" << std::string(75,'-') << std::endl; */
  /*         out << "#" << "\t Solution: " << id << std::endl; */
  /*         out << "#" << std::string(75,'-') << std::endl; */
  /*       } */

  /*       for(unsigned int c=0; c<this->_totalNCoeff; c++) */
  /*       { */
  /*         unsigned int myC = (c-coeffStart)/commSize; */
  /*         std::vector<double> myCoeff(_solSize[id],0.); */
  /*         libMesh::Parallel::MessageTag tag(c); */
  /*         libMesh::Parallel::Status stat; */

  /*         if (c%commSize == myRank) */
  /*         { */
  /*           if (myRank == 0) */
  /*           { */
  /*             for(unsigned int comp=0; comp < (_coefficients[id])[myC].size(); comp++) */
  /*               out << std::setprecision(5) << std::scientific */ 
  /*                 << (_coefficients[id])[myC](comp) << " " ; */
  /*             out << std::endl; */
  /*           } */
  /*           else */
  /*           { */
  /*             for(unsigned int comp=0; comp<myCoeff.size(); comp++) */
  /*               myCoeff[comp]= _coefficients[id][myC](comp) ; */ 
  /*             _comm.send(0,myCoeff,tag); */

  /*           } */
  /*         } */
  /*         else */
  /*         { */
  /*           if (myRank == 0) */
  /*           { */
  /*             stat=_comm.receive( c%commSize, myCoeff, tag); */
  /*             for(unsigned int comp=0; comp < myCoeff.size(); comp++) */
  /*               out << std::setprecision(5) << std::scientific */ 
  /*                 << myCoeff[comp] << " " ; */
  /*             out << std::endl; */
  /*           } */
  /*         } */


  /*       } */
  /*     } */


  /*   } */

/********************************************//**
 * \brief print coefficient values
 ***********************************************/
  /* template<class T_S> */
  /*   void SurrogateModel<T_S>::printCoefficients( */ 
  /*       std::string solutionName, */
  /*       std::ostream& out ) */ 
  /*   { */
  /*     std::vector<std::string> solutionsToPrint(1,solutionName); */
  /*     printCoefficients(solutionsToPrint, out); */
  /*   } */

/********************************************//**
 * \brief print coefficient values
 ***********************************************/
  /* template<class T_S> */
  /*   void SurrogateModel<T_S>::printCoefficients( */ 
  /*       std::ostream& out ) */ 
  /*   { */
  /*     std::vector< std::string > solutionsToPrint ; */
      
  /*     typename std::map< std::string, std::vector<NumericVector<Number> > >::iterator id; */
  /*     for (id=_coefficients.begin(); id!=_coefficients.end(); id++) */
  /*       solutionsToPrint.push_back( id->first ) ; */
  /*     printCoefficients(solutionsToPrint, out); */
  /*   } */


/********************************************//**
 * \brief basic evaluation routine, based on more general (virtual) evaluate
 * routine
 ***********************************************/
  template<class T_S> 
    std::vector<Number>  SurrogateModel<T_S>::evaluate( 
        std::string solutionName,  ///< solution to return
        T_S& parameterValues     ///< parameter values to evaluate*/
        ) 
    {
      std::vector< std::string > solutionsToGet(1,solutionName);
      std::map< std::string, std::vector<Number>  > solutionVectors
        = this->evaluate( solutionsToGet, parameterValues ) ;

      return solutionVectors[solutionName];
    }

/********************************************//**
 * \brief basic evaluation routine, based on more general (virtual) evaluate
 * routine
 ***********************************************/
  template<class T_S> 
    std::map<std::string, std::vector<Number> > SurrogateModel<T_S>::evaluate( 
        T_S& parameterValues     ///< parameter values to evaluate*/
        ) 
    {
      std::vector< std::string > solutionsToGet ;
      
      std::set<std::string>::iterator id;
      for (id=_solutionNames.begin(); id!=_solutionNames.end(); id++)
        solutionsToGet.push_back( *id ) ;

      return evaluate( solutionsToGet, parameterValues ) ;    
    }

/********************************************//**
 * \brief calculate l2norm for single function
 * calls virtual l2Norm( solutionNames) routine
 ***********************************************/
  /* template<class T_S> */ 
  /*   NumericVector<Number>  SurrogateModel<T_S>::l2Norm( */ 
  /*       std::string solutionName  ///< solution to return */
  /*       ) */ 
  /*   { */
  /*     std::vector< std::string > solutionsToGet(1,solutionName); */
  /*     std::map< std::string, NumericVector<Number>  > solutionVectors */
  /*       = this->l2Norm( solutionsToGet ) ; */

  /*     return solutionVectors[solutionName]; */
  /*   } */

/********************************************//**
 * \brief calculate l2norm for all functions
 * calls virtual l2Norm( solutionNames) routine
 ***********************************************/
  /* template<class T_S> */ 
  /*   std::map<std::string, NumericVector<Number> > SurrogateModel<T_S>::l2Norm( ) */ 
  /*   { */
  /*     std::vector< std::string > solutionsToGet ; */
      
  /*     typename std::map< std::string, PhysicsFunction<T_S>* >::iterator id; */
  /*     for (id=_solutionFunction.begin(); id!=_solutionFunction.end(); id++) */
  /*       solutionsToGet.push_back( id->first ) ; */

  /*     return l2Norm( solutionsToGet) ; */    
  /*   } */
  
  /********************************************//**
   * \brief get mean of surrogate model
   *      returns first coeff vector
   ***********************************************/
  /* template< class T_S> */
  /*    std::map< std::string, NumericVector<Number>  > SurrogateModel<T_S>::mean( ) */
  /*   { */
  /*     std::map< std::string, NumericVector<Number>  > meanCoefficients; */

  /*     // some initialization */
  /*     unsigned int myRank = _comm.rank(); */
  /*     unsigned int commSize = _comm.size(); */
  /*     unsigned int coeffStart = std::min(_comm.rank(), _totalNCoeff-1); */

  /*     // mean is the first coefficient */
  /*     unsigned int coeff = 0 ; */

  /*     // some initialization */
  /*     unsigned int myC = (coeff-coeffStart)/commSize; */
  /*     libMesh::Parallel::MessageTag tag(coeff); */
  /*     libMesh::Parallel::Status stat; */

  /*     std::set<std::string>::iterator id = _solutionNames.begin(); */
  /*     for (; id!=_solutionNames.end(); ++id) */
  /*     { */
  /*       std::vector<double> myCoeff(_solSize[*id],0.) ; */
  /*       NumericVector<Number>  solutionCoeff( _solSize[*id] ); */

  /*       if (coeff%commSize == myRank) */
  /*       { */
  /*         // set value of mean coefficient */
  /*         for(unsigned int comp=0; comp<myCoeff.size(); comp++) */
  /*           myCoeff[comp]= _coefficients[*id][myC](comp) ; */ 
  /*       } */

  /*       // if I don't have coeff 0 receive it if I do send it */
  /*       _comm.broadcast(myCoeff, coeff%commSize); */

  /*       // insert into global set */
  /*       for(unsigned int i=0; i<solutionCoeff.size(); i++) */
  /*         solutionCoeff(i) = myCoeff[i]; */

  /*       meanCoefficients.insert( */
  /*           std::pair< std::string, NumericVector<Number> >( *id, solutionCoeff ) */
  /*           ); */
  /*     } */

  /*     return meanCoefficients; */
  /*   } */
  

  /********************************************//**
   * \brief 
   ***********************************************/
  /* template< class T_S> */
  /*   double SurrogateModel<T_S>::l2NormDifference( */
  /*       SurrogateModel<T_S>& comparisonModel, */
  /*       std::string solutionName ) */
  /*   { */
  /*     std::cerr << "\n\t ERROR: l2NormDifference( ... )" */
  /*       << " is not implemented in derived class\n" */ 
  /*       << std::endl; */
  /*     exit(1); */
  /*   } */


}
#endif //SURROGATE_MODEL_H
