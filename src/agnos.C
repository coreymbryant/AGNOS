
#include <iostream>
#include <cstring>
#include <GetPot>

void run (GetPot& input);

int main(int argc, char* argv[])
{

  if( argc < 2 )
  {
    std::cerr << "\n\t ERROR: must specify input file.\n" << std::endl;
    exit(1);
  }


  GetPot input( argv[1] );

  run( input );

  return 0;
}


