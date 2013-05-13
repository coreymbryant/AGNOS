
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

      virtual void run( ) = 0;

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
    return;
  }

}


#endif // DRIVER_H


