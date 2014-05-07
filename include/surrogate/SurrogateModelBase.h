
#ifndef SURROGATE_MODEL_BASE_H
#define SURROGATE_MODEL_BASE_H

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
   * \brief Abstract SurrogateModelBase class.
   *
   * Used to derive a physics based surrogate model (SurrogateModel) and a
   * simple evaluator surrogate model (SurrogateEvaluator). 
   ***********************************************/
  template<class T_S,class T_P>
  class SurrogateModelBase
  {
    public: 

      /** Constructor */
      SurrogateModelBase(
          const Communicator&               comm,
          const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
          const std::vector<unsigned int>&          order,
          std::set<std::string> computeSolutions 
          );
      /** Destructor */
      virtual ~SurrogateModelBase();

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
          SurrogateModelBase<T_S,T_P>& comparisonModel,
          std::string solutionName
          ) ;

      /** set parameters object  */
      void setParameters( 
          std::vector<std::shared_ptr<AGNOS::Parameter> >& parameters );
      /** return reference to parameters object */
      const std::vector<std::shared_ptr<AGNOS::Parameter> >   
        getParameters( ) const { return _parameters; }

      /** reference to locally stored coefficients */
      /* const std::map< std::string, LocalMatrix>   getLocalCoefficients() const */
      /* { return _coefficients; } */
      /** reference to all coefficients */
      const std::map< std::string, LocalMatrix>   getCoefficients() ;

      /** print the coefficient vectors */
      void printCoefficients( 
        std::vector<std::string> solutionNames,
        std::ostream& out ) ;
      /** print the coefficient vectors */
      void printCoefficients( std::string solutionName, std::ostream& out ) ;
      /** print the coefficient vectors */
      void printCoefficients( std::ostream& out ) ;

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

      /** Parameter dimension */
      unsigned int dimension(){ return _dimension; }

      /** solution names this surrogateModel is built for */
      std::set<std::string> getSolutionNames( ) const
      { return _solutionNames; }

      /** reference to sol vector sizes */
      std::map<std::string, unsigned int> getSolSize( ) const
      { return _solSize; }

      /** reference to physicsGroup*/
      const int physicsGroup() const
      { return _physicsGroup; }
      /** reference to number of physics groups*/
      const int nPhysicsGroups() const
      { return _nPhysicsGroups; }
      /** reference to groupRank */
      const int groupRank() const
      { return _groupRank; }

    protected: 

      /** reference to communicator */
      const Communicator& _comm;
      /** reference to physics group number */
      int _physicsGroup;
      /** reference to number of physics groups */
      int _nPhysicsGroups;
      /** reference to groupRank (needed to deterimine if this is a master or
       * slave node) */
      int _groupRank;

      /** expansion order */
      std::vector<unsigned int> _order;  

      /** coefficients vectors */
      std::map< std::string, std::shared_ptr<DistMatrix> >  _coefficients;

      /** total number of coefficients on all vectors */
      unsigned int _totalNCoeff;
      /** indicies of coefficients corresponding to local process */
      std::vector<unsigned int> _coeffIndices;

      /** dimension of coefficient vectors for each solution name */
      std::map< std::string, unsigned int> _solSize;

      /** reference to Parameter object */
      std::vector<std::shared_ptr<AGNOS::Parameter> > _parameters;

      /** Parameter dimension */
      unsigned int _dimension;

      /** Set of solution names that the surrogate model computes */
      std::set<std::string>                               _solutionNames;
  } ;


} // namespace AGNOS
#endif // SURROGATE_MODEL_BASE_H
