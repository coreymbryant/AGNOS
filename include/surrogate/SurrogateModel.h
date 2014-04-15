
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
          const std::vector<unsigned int>&  order,
          std::set<std::string> computeSolutions = std::set<std::string>()
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

      /** Refine the surrogate model.  */
      virtual void refine( );
      virtual void refine( const std::vector<unsigned int>& increase ) ;

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

      /** Sample the surrogate model N times and place results in sampleVec  */
      virtual void sample( 
          std::string solutionName, unsigned int N, std::vector<T_P>& sampleVec
          );

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


}
#endif //SURROGATE_MODEL_H
