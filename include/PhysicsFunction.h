
#ifndef PHYSICS_FUNCTION_H
#define PHYSICS_FUNCTION_H

#include "PhysicsModel.h"

namespace AGNOS
{

  template<class T_S, class T_P>
  class PhysicsFunction
  {

    public:

      PhysicsFunction(
          PhysicsModel<T_S,T_P>& physics
          );
      virtual ~PhysicsFunction();

      virtual void compute( 
          const T_S& paramVector,
          T_P& imageVector
          ) = 0;

    protected:
      PhysicsModel<T_S,T_P>& m_physics;

  };

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunction<T_S,T_P>::PhysicsFunction( 
        PhysicsModel<T_S,T_P>& physics
        )
    : m_physics(physics)
    {
    }

/********************************************//**
 * \brief 
 *
 * 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunction<T_S,T_P>::~PhysicsFunction( )
    {
    }

  

}


#endif // PHYSICS_FUNCTION_H
