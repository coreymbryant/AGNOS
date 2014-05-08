
#ifndef SURROGATE_EVALUATOR_H
#define SURROGATE_EVALUATOR_H

#include "SurrogateModelBase.h"

namespace AGNOS
{

  /********************************************//**
   * \brief Evaluator surrogate model
   *
   * Class built on base SurrogateModel class that should be used only for
   * constructing an evaluator based on provided coefficients, etc. Ideal for
   * external use of a final surrogate model.
   ***********************************************/
  template<class T_S,class T_P>
  class SurrogateEvaluator : virtual public SurrogateModelBase<T_S,T_P> 
  {
    public: 
      /** Constructor */
      SurrogateEvaluator(
          const Communicator&               comm,
          const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
          const std::vector<unsigned int>&          order,
          std::set<std::string> computeSolutions 
          ) ;

      /** Destructor */
      ~SurrogateEvaluator( );

      /** build the surrogate model */
      virtual void build( 
          std::vector< std::vector<unsigned int> > indexSet,
          std::map< std::string, std::shared_ptr<DistMatrix> >  coefficients
          ) ; 

      /********************************************//**
       * \brief build surrogate model
       *
       * must be redefined in derived classes because type of surrogate impacts
       * how coefficients are stored accross multiple processors
       ***********************************************/
      virtual void build( 
          std::vector< std::vector<unsigned int> > indexSet,
          std::map< std::string, LocalMatrix >  coefficients
          ) = 0; 


    protected:

  };
  

} // namespace AGNOS
#endif // SURROGATE_EVALUATOR_H
