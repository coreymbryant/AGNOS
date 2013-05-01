
#ifndef DRIVER_H
#define DRIVER_H


#include <iostream>
#include <stdio.h>
#include <assert.h>
#include "TensorProductPseudoSpectral.h"
#include "CatenaryModel.h"

namespace AGNOS
{
  typedef std::vector<double> T_S ;
  typedef std::vector<double> T_P ;
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

  Driver::Driver( )
    : m_surrogateOrder(1)
  {
  }

  Driver::Driver(int surrogateOrder )
    : m_surrogateOrder(surrogateOrder)
  {
  }

  Driver::~Driver( )
  {
  }

  // an initial driver run routine for testing
  void Driver::run( )
  {

    unsigned int dimension = 2;

    std::vector<Parameter*> myParameters(dimension, new Parameter()); 
    myParameters[0] = new Parameter(UNIFORM,0.0,1.0 );
    myParameters[1] = new Parameter(UNIFORM,0.0,1.0 );

    std::vector<unsigned int> myOrder(dimension,1);

    PhysicsModel<T_S,T_P>* myPhysics = new CatenaryModel<T_S,T_P>( ) ;

    typedef T_P* (PhysicsModel<T_S,T_P>::*PhysicsMemberFcn)(T_S&);
    PhysicsMemberFcn physicsFunction = &PhysicsModel<T_S,T_P>::solvePrimal;
    /* PhysicsFunc myFunc = myPhysics->solvePrimal ; */


    TensorProductPseudoSpectral<T_S,T_P>* mySurrogate = new 
      TensorProductPseudoSpectral<T_S,T_P>(
          physicsFunction, // TODO THIS STILL WONT WORK
          myParameters, 
          myOrder  
          );


    mySurrogate->getQuadRule()->printQuadPoints( );
    mySurrogate->getQuadRule()->printQuadWeights( );

    return;
  }

}


#endif // DRIVER_H


