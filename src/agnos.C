
#include <mpi.h>
#include "agnosDefines.h"
#include "Driver.h"

int main(int argc, char* argv[])
{

  if( argc < 2 )
  {
    std::cerr << "\n\t ERROR: must specify AGNOS input file.\n" << std::endl;
    exit(1);
  }

  std::string file_name = argv[1];

  GetPot inputfile( file_name);


  MPI_Init(&argc,&argv);
  Communicator comm(MPI_COMM_WORLD);

  {
    MPI_Comm myComm;
    int mpiSplit =  
      MPI_Comm_split( MPI_COMM_WORLD, comm.rank(), 0, &myComm);

    LibMeshInit libmesh_init(argc, argv, myComm);
    
    libMesh::Parallel::Communicator physicsComm(myComm);

    AGNOS::Driver agnos( comm, physicsComm, inputfile );

    agnos.run( );
    
    /* LibMeshInit libmesh_init(argc, argv); */
  }

#ifndef AGNOS_USING_MPICH
  int ierr;
  comm.barrier();
  MPI_Finalize();
#endif
  return 0;
}

