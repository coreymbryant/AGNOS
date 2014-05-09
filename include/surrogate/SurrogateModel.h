
#ifndef SURROGATE_MODEL_H 
#define SURROGATE_MODEL_H 

#include "SurrogateModelBase.h"

namespace AGNOS
{

  enum SurrogateModelType{
    PSEUDO_SPECTRAL_TENSOR_PRODUCT=0,
    PSEUDO_SPECTRAL_SPARSE_GRID,
    PSEUDO_SPECTRAL_MONTE_CARLO,
    COLLOCATION };


  /********************************************//**
   * \brief Surrogate model class
   *
   * Allows for derivation of Collocation and Pseudospectral surrogate models.
   * PhysicsModel must be provided
   ***********************************************/
  template<class T_S, class T_P>
  class SurrogateModel : virtual public SurrogateModelBase<T_S,T_P>
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
          std::shared_ptr<SurrogateModelBase<T_S,T_P> > primarySurrogate, 
          std::vector<unsigned int> increaseOrder = std::vector<std::string>(),
          unsigned int multiplyOrder = 1,
          std::set<std::string> evaluateSolutions = std::vector<std::string>(),
          std::set<std::string> computeSolutions = std::vector<std::string>()
          );

      /** default destructor */
      virtual ~SurrogateModel( ); 

      /** Initialization routine */
      virtual void initialize( ) = 0 ;

      virtual void refineUniformly( ) ;
      virtual void refine( const std::vector<unsigned int>& increase ) ;

    protected: 
      
      /** expansion order */
      std::vector<unsigned int> _increaseOrder ;
      unsigned int              _multiplyOrder ;

      /** Set of solution names that the surrogate model usese primary surrogate
       * for. Only meaningful in secondary constructor case */
      std::set<std::string>                               _evalNames;

      /** Primary surrogate to use in evaluation for secondary surrogate
       * construciton */
      std::shared_ptr<SurrogateModelBase<T_S,T_P> >      _evalSurrogate ;

      /** Data structure to hold evalSurrogate evaluations, to be used in
       * surrogate construction */
      std::vector< std::map< std::string,T_P> > _primaryEvaluations;
      

  }; //SurrogateModel class


}
#endif //SURROGATE_MODEL_H
