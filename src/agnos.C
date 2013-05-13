
#include <vector>
#include <iostream>
#include <cstring>
#include <GetPot>
#include "Driver.h"

int main(int argc, char* argv[])
{

  if( argc < 2 )
  {
    std::cerr << "\n\t ERROR: must specify input file.\n" << std::endl;
    exit(1);
  }


  GetPot input( argv[1] );


  AGNOS::Driver agnos( input ) ;

  agnos.run();

  return 0;
}


