
#ifndef DRIVER_H
#define DRIVER_H


#include <iostream>
#include <stdio.h>
#include <assert.h>
// TODO find a better way to include all necessary submodels?
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

  /********************************************//**
   * \brief Base class for driving the construction of SurrogateModel objects
   *
   * This Base class is designed to be a non-adaptive (i.e. single construction)
   * method for surrogate construction.  Several different adaptive strategies
   * can be implemented as derived classes of this routine. 
   *
   * TODO this should be defined based on user inputs. Or at the very least
   * require a user provided derived class to be defined for the particular
   * application.
   * 
   *
   * 
   ***********************************************/
  class Driver
  { 

    public:

      Driver( );
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

    std::vector<Parameter*> myParameters(
        dimension, 
        new Parameter(UNIFORM, 1.0,3.0)
        ); 

    std::vector<unsigned int> myOrder(dimension,1);
    myOrder.front() = 1;

    PhysicsCatenary<T_S,T_P>* myPhysics = new PhysicsCatenary<T_S,T_P>( ) ;

    PhysicsFunction<T_S,T_P>* myPhysicsFunction =
      new PhysicsFunctionPrimal<T_S,T_P>( *myPhysics ) ;

    PseudoSpectralTensorProduct<T_S,T_P>* mySurrogate = new 
      PseudoSpectralTensorProduct<T_S,T_P>(
          *myPhysicsFunction, 
          myParameters, 
          myOrder  
          );

    mySurrogate->build( );

    std::vector<T_P> myCoeff = mySurrogate->getCoefficients( );
    for(unsigned int coeff=0; coeff<myCoeff.size(); coeff++)
        std::cout << std::setprecision(5) << std::scientific 
          << "coeff[" << coeff << "](0) = " 
          << myCoeff[coeff](0) << std::endl;


    T_S paramValue(dimension);
    paramValue(0) = 1.5;

    T_P testValue = mySurrogate->evaluate( paramValue ) ;

    std::cout << "testValue(0) = " << testValue(0) << std::endl;

    return;
  }

}


#endif // DRIVER_H

