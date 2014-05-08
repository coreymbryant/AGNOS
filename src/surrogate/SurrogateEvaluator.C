
#include "SurrogateEvaluator.h"

namespace AGNOS
{
  /********************************************//**
   * \brief Constructor
   ***********************************************/
  template<class T_S,class T_P>
    SurrogateEvaluator<T_S,T_P>::SurrogateEvaluator(
        const Communicator&               comm,
        const std::vector<std::shared_ptr<AGNOS::Parameter> >&     parameters,
        const std::vector<unsigned int>&          order,
        std::set<std::string> computeSolutions 
        )
    : SurrogateModelBase<T_S,T_P>(comm,parameters,order,computeSolutions)
    {
    }

  /********************************************//**
   * \brief Destructor
   ***********************************************/
  template<class T_S,class T_P>
    SurrogateEvaluator<T_S,T_P>::~SurrogateEvaluator()
    { 
      this->_coefficients.clear(); 
      this->_indexSet.clear();
    }



  /********************************************//**
   * \brief Build routine. Sets indexSet and Coefficients 
   ***********************************************/
  template<class T_S,class T_P>
    void SurrogateEvaluator<T_S,T_P>::build( 
        std::vector< std::vector<unsigned int> > indexSet,
        std::map< std::string, std::shared_ptr<DistMatrix> >  coefficients
        ) 
    {
      this->_coefficients = coefficients;
      this->_indexSet = indexSet;

      return;
    }


  template class 
    SurrogateEvaluator<libMesh::DenseVector<double>, libMesh::DenseVector<double> >;
} //namespace AGNOS
