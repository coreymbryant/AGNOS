
#include "agnosDefines.h"

void run ( const libMesh::Parallel::Communicator& comm, GetPot& input);

int main(int argc, char* argv[])
{

  if( argc < 2 )
  {
    std::cerr << "\n\t ERROR: must specify AGNOS input file.\n" << std::endl;
    exit(1);
  }


  MPI_Init(&argc,&argv);
  const libMesh::Parallel::Communicator comm(MPI_COMM_WORLD);


  GetPot inputfile( argv[1] );

  MPI_Comm myComm;
  int mpiSplit =  
    MPI_Comm_split( MPI_COMM_WORLD, comm.rank(), 0, &myComm);

  LibMeshInit libmesh_init(argc, argv, myComm);
  

  run( comm,  inputfile );
  


  /* MPI_Finalize(); */
  return 0;
}


