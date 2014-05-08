
#ifndef EVALUATOR_PSEUDO_SPECTRAL_H
#define EVALUATOR_PSEUDO_SPECTRAL_H

#include "SurrogateEvaluator.h"

namespace AGNOS
{

  /********************************************//**
   * \brief Tensor product evaluator. 
   *
   * An evaluator class based on data from a PseudoSpectralPseudoSpectral
   * surrogate model.
   * 
   ***********************************************/
  template<class T_S,class T_P>
  class EvaluatorPseudoSpectral : public SurrogateEvaluator<T_S,T_P>
  {
    public: 
      /** Constructor */
      EvaluatorPseudoSpectral(
          const Communicator&               comm,
          const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
          const std::vector<unsigned int>&          order,
          std::set<std::string> computeSolutions 
          ) ;

      /** Destructor */
      ~EvaluatorPseudoSpectral( );

      /** build the surrogate model */
      using SurrogateEvaluator<T_S,T_P>::build;
      virtual void build( 
          std::vector< std::vector<unsigned int> > indexSet,
          std::map< std::string, LocalMatrix >  coefficients
          ) ; 

      /** Evaluate basis polys */
      std::vector<double>       evaluateBasis( 
          const std::vector< std::vector<unsigned int> >& indexSet,
          T_S& parameterValues ) const;
      /** Evaluation routine */
      virtual std::map<std::string, T_P> evaluate( 
          std::set< std::string >  solutionNames,
          T_S&                        parameterValues,
          bool saveLocal = true /**< save solution locally after evaluation*/
          ) const ;


    protected:

  }; // class EvaluatorPseudoSpectral

} // namespace AGNOS
#endif // EVALUATOR_PSEUDO_SPECTRAL_H
