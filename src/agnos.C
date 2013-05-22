
#include "agnosDefines.h"

void run ( const libMesh::Parallel::Communicator& comm, GetPot& input);

int main(int argc, char* argv[])
{

  if( argc < 2 )
  {
    std::cerr << "\n\t ERROR: must specify input file.\n" << std::endl;
    exit(1);
  }

  MPI_Init(&argc,&argv);
  
  const libMesh::Parallel::Communicator comm(MPI_COMM_WORLD);

  GetPot input( argv[1] );

  run( comm,  input );

  MPI_Finalize();
  return 0;
}


