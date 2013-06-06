
#include "agnosDefines.h"

void run ( const libMesh::Parallel::Communicator& comm, GetPot& input);

int main(int argc, char* argv[])
{

  if( argc < 2 )
  {
    std::cerr << "\n\t ERROR: must specify AGNOS input file.\n" << std::endl;
    exit(1);
  }


  /* MPI_Init(&argc,&argv); */
  /* const libMesh::Parallel::Communicator comm(MPI_COMM_WORLD); */

  GetPot inputfile( argv[1] );

  LibMeshInit libmesh_init(argc, argv);

  run( libMesh::Parallel::Communicator_World,  inputfile );
  


  /* MPI_Finalize(); */
  return 0;
}


