
#ifndef PHYSICS_FUNCTION_TOTAL_ERRORH
#define PHYSICS_FUNCTION_TOTAL_ERRORH

#include <map>

#include "PhysicsModel.h"
#include "SurrogateModel.h"

namespace AGNOS
{


/********************************************//**
 * \brief Total Error function based on a PhysicsModel object
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunctionTotalError : public PhysicsFunction<T_S,T_P>
  {
    public:
      PhysicsFunctionTotalError( PhysicsModel<T_S,T_P>* physics,
          SurrogateModel<T_S,T_P>* solutionSurrogate );
      void compute( const T_S& paramVector, T_P& imageVector) ;

      void computeData( std::vector<T_S> integrationPoints );
      void setData( unsigned int currentIndex );

    protected:
      SurrogateModel<T_S,T_P>* m_solutionSurrogate;
      std::vector<T_P> m_primalData;
      std::vector<T_P> m_adjointData;

  };
/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunctionTotalError<T_S,T_P>::PhysicsFunctionTotalError( 
        PhysicsModel<T_S,T_P>* physics,
        SurrogateModel<T_S,T_P>* solutionSurrogate 
        ) : m_solutionSurrogate(solutionSurrogate)
    { 
      this->m_name = "error"; 
      this->m_physics = physics;
    }

/********************************************//**
 * \brief compute primal and adjoint for all integration points. Will be needed
 * in the compute(...) routine
 ***********************************************/
  template<class T_S, class T_P> 
    void PhysicsFunctionTotalError<T_S,T_P>::computeData( 
        std::vector<T_S> integrationPoints )
    {
      m_primalData.clear();
      m_adjointData.clear();
      
      // make sure solution and adjoint are present in provided surrogate model
      if (m_solutionSurrogate->getSolutionNames().count("primal") == 0)
      {
        std::cerr << "\n\t ERROR: primal solution not present in "
          << " surrogate model needed to compute TotalError. \n "
          << std::endl;
        exit(1);
      }
      if (m_solutionSurrogate->getSolutionNames().count("adjoint") == 0)
      {
        std::cerr << "\n\t ERROR: adjoint solution not present in "
          << " surrogate model needed to compute TotalError. \n "
          << std::endl;
        exit(1);
      }

      // set solution names
      std::vector<std::string> names;
      names.push_back("adjoint");
      names.push_back("primal");
      
      // generate data for solution and adjoint
      for (unsigned int i=0; i<integrationPoints.size(); i++)
      {
        std::map<std::string, T_P> surrogateEvaluations 
          = m_solutionSurrogate->evaluate(names, integrationPoints[i] );

        m_primalData.push_back( surrogateEvaluations["primal"] );
        m_adjointData.push_back( surrogateEvaluations["adjoint"] );
      }
      
      
      /* std::cout << "test: computeData(...) " << m_primalData.size() << std::endl; */
      return;
    }

/********************************************//**
 * \brief set primal and adjoint to current integration point data
 ***********************************************/
  template<class T_S, class T_P> 
    void PhysicsFunctionTotalError<T_S,T_P>::setData(
        unsigned int currentIndex )
    {
      this->m_physics->setPrimalSolution( m_primalData[currentIndex] );
      this->m_physics->setAdjointSolution( m_adjointData[currentIndex] );

      /* std::cout << "m_primalData.size(): " << m_primalData[currentIndex].size() */
      /*                                    << std::endl; */
      /* for(unsigned int i=0; i<m_primalData[currentIndex].size(); i++) */
      /*   std::cout << "m_primalData(" << i << ") = " << */
      /*     m_primalData[currentIndex](i) << std::endl; */

      /* std::cout << "primal.size(): " << this->m_physics->getPrimalSolution()->size() */
      /*                                    << std::endl; */
      /* for(unsigned int i=0; i<this->m_physics->getPrimalSolution()->size(); i++) */
      /*   std::cout << "primal(" << i << ") = " << */
      /*     (*this->m_physics->getPrimalSolution())(i) << std::endl; */

      /* std::cout << "m_adjointData.size(): " << m_adjointData[currentIndex].size() */
      /*                                    << std::endl; */
      /* for(unsigned int i=0; i<m_adjointData[currentIndex].size(); i++) */
      /*   std::cout << "m_adjointData(" << i << ") = " << */
      /*     m_adjointData[currentIndex](i) << std::endl; */

      /* std::cout << "adjoint.size(): " << this->m_physics->getAdjointSolution()->size() */
      /*                                    << std::endl; */
      /* for(unsigned int i=0; i<this->m_physics->getAdjointSolution()->size(); i++) */
      /*   std::cout << "adjoint(" << i << ") = " << */
      /*     (*this->m_physics->getAdjointSolution())(i) << std::endl; */


      /* std::cout << "test: setData(...) " << std::endl; */
      return;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    void PhysicsFunctionTotalError<T_S,T_P>::compute(
        const T_S& paramVector,
        T_P& imageVector
        )
    {

      /* std::cout << "test: compute beginning " << std::endl; */

      /* std::cout << "primal.size(): " << this->m_physics->getPrimalSolution()->size() */
      /*                                    << std::endl; */
      /* for(unsigned int i=0; i<this->m_physics->getPrimalSolution()->size(); i++) */
      /*   std::cout << "primal(" << i << ") = " << */
      /*     (*this->m_physics->getPrimalSolution())(i) << std::endl; */
      /* std::cout << "adjoint.size(): " << this->m_physics->getAdjointSolution()->size() */
      /*                                    << std::endl; */
      /* for(unsigned int i=0; i<this->m_physics->getAdjointSolution()->size(); i++) */
      /*   std::cout << "adjoint(" << i << ") = " << */
      /*     (*this->m_physics->getAdjointSolution())(i) << std::endl; */

      T_P tempPrimal(*(this->m_physics->getPrimalSolution()) ) ;
      T_P tempAdjoint(*(this->m_physics->getAdjointSolution()) ) ;

      /* std::cout << "test: compute post temp vectors " << std::endl; */
      imageVector = this->m_physics->estimateError(
          paramVector, tempPrimal, tempAdjoint);

      /* std::cout << "test: compute end " << std::endl; */
      return;
    }


}


#endif // PHYSICS_FUNCTION_TOTAL_ERRORH
