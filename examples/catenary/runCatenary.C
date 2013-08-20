
#include "agnosDefines.h"
#include "Driver.h"

#include "PhysicsCatenary.h"

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

  AGNOS::Driver myDriver(comm,comm,input);

  myDriver.run( );


  return;
}


