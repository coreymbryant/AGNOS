
#ifndef PHYSICS_FUNCTION_H
#define PHYSICS_FUNCTION_H

#include "PhysicsModel.h"
#include <assert.h>

namespace AGNOS
{

/********************************************//**
 * \brief Base physics function class
 *
 * This class provides a means for accessing primal and adjoint solutions from
 * the PhysicsModel class. It also allows the definition of simple
 * (user-defined) functions that are independent of a PhysicsModel object. 
 *
 * It provides a convenient, uniform, interface for a SurrogateModel object to
 * call, namely the compute() function. 
 * 
 * 
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunction
  {

    public:

      PhysicsFunction();
      virtual ~PhysicsFunction();

      virtual void compute( 
          const T_S& paramVector,
          T_P& imageVector
          ) = 0;

      unsigned int getImageSize() const ;

    protected:
      unsigned int m_imageSize;

  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunction<T_S,T_P>::PhysicsFunction( ) 
    { }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunction<T_S,T_P>::~PhysicsFunction( )
    { }
  
/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    unsigned int PhysicsFunction<T_S,T_P>::getImageSize() const 
    {
      return m_imageSize;
    }

/********************************************//**
 * \brief Basic function independent of a PhysicsModel object
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunctionSimple : public PhysicsFunction<T_S,T_P>
  {
    public:
      PhysicsFunctionSimple( 
          T_P (*myFunction)(const T_S&),
          unsigned int imageSize 
        )
        : m_myFunction(myFunction)
    {
      this->m_imageSize(imageSize) ;
    }

      void compute( const T_S& paramVector, T_P& imageVector)
      {
        if (imageVector.size() != this->m_imageSize )
        {
          std::cout << "\n ERROR: returned physics vector does not match"
            << " predefined imageVector size\n \n " ;
        assert(0);
        }

        imageVector = m_myFunction(paramVector) ;
        return ;
      };

    protected:
      T_P (*m_myFunction)(const T_S&);
  };


/********************************************//**
 * \brief Primal solution function based on a PhysicsModel object
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunctionPrimal : public PhysicsFunction<T_S,T_P>
  {
    public:
      PhysicsFunctionPrimal( 
          PhysicsModel<T_S,T_P>& physics,
          unsigned int imageSize
          );
      void compute( const T_S& paramVector, T_P& imageVector) ;
    protected:
      PhysicsModel<T_S,T_P>& m_physics;
  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunctionPrimal<T_S,T_P>::PhysicsFunctionPrimal( 
        PhysicsModel<T_S,T_P>& physics,
        unsigned int imageSize
        )
        : m_physics(physics)
    {
      this->m_imageSize = imageSize ;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    void PhysicsFunctionPrimal<T_S,T_P>::compute(
        const T_S& paramVector,
        T_P& imageVector
        )
    {
      if (imageVector.size() != this->m_imageSize )
      {
        std::cout << "\n ERROR: returned physics vector does not match"
          << " predefined imageVector size\n \n " ;
        assert(0);
      }

      imageVector = this->m_physics.solvePrimal(paramVector);
      this->m_physics.setPrimalSolution(imageVector);
      return ;
    }
  

/********************************************//**
 * \brief Adjoint solution function based on a PhysicsModel object
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunctionAdjoint : public PhysicsFunction<T_S,T_P>
  {
    public:
      PhysicsFunctionAdjoint( 
          PhysicsModel<T_S,T_P>& physics,
          unsigned int imageSize
          );
      void compute( const T_S& paramVector, T_P& imageVector) ;
    protected:
      PhysicsModel<T_S,T_P>& m_physics;
  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunctionAdjoint<T_S,T_P>::PhysicsFunctionAdjoint( 
        PhysicsModel<T_S,T_P>& physics,
        unsigned int imageSize
        )
        : m_physics(physics)
    {
      this->m_imageSize = imageSize ;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    void PhysicsFunctionAdjoint<T_S,T_P>::compute(
        const T_S& paramVector,
        T_P& imageVector
        )
    {
      if (imageVector.size() != this->m_imageSize )
      {
        std::cout << "\n ERROR: returned physics vector does not match"
          << " predefined imageVector size\n \n " ;
        assert(0);
      }
      imageVector = this->m_physics.solveAdjoint(paramVector);
    }
  


/********************************************//**
 * \brief QoI solution function based on a PhysicsModel object
 ***********************************************/
  template<class T_S, class T_P>
  class PhysicsFunctionQoi : public PhysicsFunction<T_S,T_P>
  {
    public:
      PhysicsFunctionQoi( 
          PhysicsModel<T_S,T_P>& physics,
          unsigned int imageSize
          );
      void compute( const T_S& paramVector, T_P& imageVector) ;
    protected:
      PhysicsModel<T_S,T_P>& m_physics;
  };

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    PhysicsFunctionQoi<T_S,T_P>::PhysicsFunctionQoi( 
        PhysicsModel<T_S,T_P>& physics,
        unsigned int imageSize
        )
        : m_physics(physics)
    {
      this->m_imageSize = imageSize ;
    }

/********************************************//**
 * \brief 
 ***********************************************/
  template<class T_S,class T_P>
    void PhysicsFunctionQoi<T_S,T_P>::compute(
        const T_S& paramVector,
        T_P& imageVector
        )
    {
      if (imageVector.size() != this->m_imageSize )
      {
        std::cout << "\n ERROR: returned physics vector does not match"
          << " predefined imageVector size\n \n " ;
        // TODO proper way to make this fail
        assert(0);
      }
      imageVector = this->m_physics.evaluateQoi(paramVector);
    }



}


#endif // PHYSICS_FUNCTION_H
