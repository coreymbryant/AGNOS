
#include "agnosDefines.h"
#include "DriverPhysics.h"

#include "PhysicsCatenaryLibmesh.h"

/********************************************//**
 * \brief User must define a specific run function to provide physics class and
 * initialize driver. Everthing else will be handled behind the scenes. 
 *
 * It is up to user how to provide physics initialization data (from input file)
 * to PhysicsModel constructor.
 *
 * The following is an example run( input ) routine for a simple scalar valued
 * PhysicsModel class defined in PhysicsCatenary.
 * 
 ***********************************************/

void run ( const Communicator& comm, GetPot& input )
{

  AGNOS::PhysicsModel<T_S,T_P>* myPhysics = 
    new AGNOS::PhysicsCatenaryLibmesh<T_S,T_P>( 
      comm, input );


  // this should be the same for any user 
  // (unless DriverFunction is used)
  AGNOS::DriverPhysics myDriver( comm, myPhysics, input );

  myDriver.run( );

  delete myPhysics;

  return;
}


