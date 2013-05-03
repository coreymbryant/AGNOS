
#ifndef DRIVER_H
#define DRIVER_H


#include <iostream>
#include <stdio.h>
#include <assert.h>
#include "PseudoSpectralTensorProduct.h"
#include "PhysicsCatenary.h"
#include "PhysicsFunction.h"

#include "libmesh/dense_vector.h"


namespace AGNOS
{
  /* typedef std::vector<double> T_P ; */
  /* typedef std::vector<double> T_S ; */
  typedef libMesh::DenseVector<double> T_P ;
  typedef libMesh::DenseVector<double> T_S ;

  class Driver
  { 

    public:

      Driver( );
      Driver(int surrogateOrder);
      virtual ~Driver( );

      void run( );

    private:
      int m_surrogateOrder;
  };

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::Driver( )
    : m_surrogateOrder(1)
  {
  }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::Driver(int surrogateOrder )
    : m_surrogateOrder(surrogateOrder)
  {
  }

/********************************************//**
 * \brief 
 * 
 ***********************************************/
  Driver::~Driver( )
  {
  }

/********************************************//**
 * \brief an initial driver run routine for testing
 * 
 ***********************************************/
  void Driver::run( )
  {

    unsigned int dimension = 3;

    std::vector<Parameter*> myParameters(
        dimension, 
        new Parameter(UNIFORM, -1.0,1.0)
        ); 

    std::vector<unsigned int> myOrder(dimension,1);
    myOrder.back() = 2;

    PhysicsCatenary<T_S,T_P>* myPhysics = new PhysicsCatenary<T_S,T_P>( ) ;

    PhysicsFunction<T_S,T_P>* myPhysicsFunction =
      new PhysicsFunctionQoi<T_S,T_P>( *myPhysics ) ;

    PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
      PseudoSpectralTensorProduct<T_S,T_P>(
          *myPhysicsFunction, 
          myParameters, 
          myOrder  
          );

    /* mySurrogate->build( ); */


    /* mySurrogate->printIntegrationPoints( ); */
    /* mySurrogate->printIntegrationWeights( ); */
    /* mySurrogate->printIndexSet( ); */

    return;
  }

}


#endif // DRIVER_H


