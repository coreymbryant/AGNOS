
#ifndef DRIVER_H
#define DRIVER_H

namespace AGNOS
{

  class Driver
  { 

    public:

      Driver( );
      ~Driver( );

      void run( );
  };


  // an initial driver run routine for testing
  Driver::run( )
  {

    int N = 1;
    PhysicsModel* physics = new CatenaryModel( );
    SurrogateModel* surrogate = new PseudoSpectral( );
    
    return;
  }

}


#endif // DRIVER_H


