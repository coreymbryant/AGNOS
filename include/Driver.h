
#ifndef DRIVER_H
#define DRIVER_H


#include <iostream>
#include <stdio.h>
#include <assert.h>
#include "TensorProductPseudoSpectral.h"
#include "CatenaryModel.h"
#include "CombinedSolutionVector.h"


namespace AGNOS
{
  typedef std::vector<double> T_P ;
  typedef std::vector<double> T_S ;
  //typedef double* T_S ;

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

    unsigned int dimension = 1;

    std::vector<Parameter*> myParameters(dimension, new Parameter()); 
    myParameters[0] = new Parameter(UNIFORM,-1.0,1.0 );
    /* myParameters[1] = new Parameter(UNIFORM,0.0,1.0 ); */

    std::vector<unsigned int> myOrder(dimension,1);

    CatenaryModel<T_S,T_P>* myPhysics = new CatenaryModel<T_S,T_P>( ) ;

    PhysicsFunction<T_S,T_P>* myPhysicsFunction =
      new CombinedSolutionVector<T_S,T_P>( *myPhysics ) ;

    TensorProductPseudoSpectral<T_S,T_P>* mySurrogate = new 
      TensorProductPseudoSpectral<T_S,T_P>(
          *myPhysicsFunction, 
          myParameters, 
          myOrder  
          );

    /* mySurrogate->build( ); */


    mySurrogate->getQuadRule()->printQuadPoints( );
    mySurrogate->getQuadRule()->printQuadWeights( );

    return;
  }

}


#endif // DRIVER_H


