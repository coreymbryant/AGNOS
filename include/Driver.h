
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
      ~Driver( );

      void run( );

    private:
      int m_N;
  };

  Driver::Driver( )
    : m_N(1)
  {
    return ;
  }

  Driver::~Driver( )
  {
    return ;
  }

  // an initial driver run routine for testing
  void Driver::run( )
  {

    PhysicsModel<T_S,T_P>* physics = new CatenaryModel<T_S,T_P>( );
    SurrogateModel<T_S,T_P>* surrogate = new PseudoSpectral<T_S,T_P>( );

    std::cout << "Hellow world \n" ;
    
    return;
  }

}


#endif // DRIVER_H


