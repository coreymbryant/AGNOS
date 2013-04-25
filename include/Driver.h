
#ifndef DRIVER_H
#define DRIVER_H


#include "PhysicsModel.h"
#include "SurrogateModel.h"
#include "CatenaryModel.h"
#include "PseudoSpectral.h"
#include <iostream>

namespace AGNOS
{
  typedef std::vector<double> T_S ;
  typedef std::vector<double> T_P ;

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
    int d=2;

    std::vector<double> mins(d,-1.0);
    std::vector<double> maxs(d,1.0);

    PhysicsModel<T_S,T_P>* physics = new CatenaryModel<T_S,T_P>( );
    PseudoSpectral<T_S,T_P>* surrogate = new
      PseudoSpectral<T_S,T_P>(2,d,mins,maxs);

    surrogate->printQuadPoints( );
    surrogate->printQuadWeights( );

    return;
  }

}


#endif // DRIVER_H


