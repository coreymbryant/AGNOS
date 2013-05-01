
#ifndef COMBINED_SOLUTION_VECTOR_H
#define COMBINED_SOLUTION_VECTOR_H
#include "PhysicsFunction.h"

namespace AGNOS
{

  template<class T_S, class T_P>
  class CombinedSolutionVector : public PhysicsFunction<T_S,T_P>
  {
    public:
      CombinedSolutionVector(
          PhysicsModel<T_S,T_P>& physics
          );
      ~CombinedSolutionVector( );

      void compute( 
          const T_S& paramVector,
          T_P& imageVector
          ) ;


    private:
  };

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
    CombinedSolutionVector<T_S,T_P>::CombinedSolutionVector( 
        PhysicsModel<T_S,T_P>& physics
        )
    : PhysicsFunction<T_S,T_P>(physics) 
    {
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
    CombinedSolutionVector<T_S,T_P>::~CombinedSolutionVector( )
    {
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
    void CombinedSolutionVector<T_S,T_P>::compute(
        const T_S& paramVector,
        T_P& imageVector
        )
    {
      T_P* primalSol = this->m_physics.getPrimalSolution( paramVector ) ;
      T_P* adjointSol = this->m_physics.getAdjointSolution( paramVector,
          this->m_physics.getPrimalSolution() ) ;
      T_P* qoiValue = this->m_physics.evaluateQoi( paramVector,
          this->m_physics.getPrimalSolution() ) ;
      unsigned int pSize = primalSol->size();
      unsigned int aSize = adjointSol->size();
      unsigned int qSize = qoiValue->size();

      imageVector.reserve( pSize + aSize + qSize );

      for (unsigned int pIt=0; pIt < pSize; pIt++)
        imageVector[pIt] = (*primalSol)[pIt];

      for (unsigned int aIt=0; aIt < aSize; aIt++)
        imageVector[ pSize + aIt] = (*adjointSol)[aIt];

      for (unsigned int qIt=0; qIt < qSize; qIt++)
        imageVector[ pSize + aSize + qIt] = (*qoiValue)[qIt];

    }

}

#endif // COMBINED_SOLUTION_VECTOR_H
