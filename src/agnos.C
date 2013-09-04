
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

  int ierr;
  ierr = MPI_Barrier(MPI_COMM_WORLD);
  ierr = MPI_Finalize();
  return 0;
}



/* int main(int argc, char* argv[]) */
/* { */
/*   if( argc < 2 ) */
/*   { */
/*     std::cerr << "\n\t ERROR: must specify AGNOS input file.\n" << std::endl; */
/*     exit(1); */
/*   } */

/*   MPI_Init(&argc,&argv); */

/*   MPI_Comm comm_physics, driverComm; */
/*   MPI_Group physicsGroup, globalGroup, driverGroup; */
/*   int rank, physicsRank, driverRank, size; */ 

/*   MPI_Comm_rank(MPI_COMM_WORLD, &rank); */
/*   MPI_Comm_size(MPI_COMM_WORLD, &size); */



/*   GetPot inputfile( argv[1] ); */

/*   int physicsNodeSize = inputfile("parallel/physicsNodeSize",1); */
/*   /1* std::cout << "physicsNodeSize:" << physicsNodeSize << std::endl; *1/ */

/*   if (size%physicsNodeSize) */
/*   { */
/*     std::cerr << "\n\t ERROR: total number of processors must be a multiple of " */
/*       << "physicsNodeSize.\n" <<std::endl; */
/*     exit(1); */
/*   } */


/*   /1* // original group *1/ */
/*   /1* MPI_Comm_group(MPI_COMM_WORLD,&globalGroup); *1/ */

/*   /1* // physics group *1/ */
/*   /1* int group = rank/physicsNodeSize; *1/ */
/*   /1* int nGroups = size/physicsNodeSize; *1/ */

/*   /1* int** groupRanks = new int*[nGroups]; *1/ */
/*   /1* for (unsigned int i=0; i<nGroups; i++) *1/ */
/*   /1*   groupRanks[i] = new int[physicsNodeSize]; *1/ */

/*   /1* for (unsigned int i=0;i<size;i++) *1/ */
/*   /1*   groupRanks[i/physicsNodeSize][i%physicsNodeSize] = i; *1/ */

/*   /1* MPI_Group_incl( *1/ */
/*   /1*     globalGroup, *1/ */ 
/*   /1*     physicsNodeSize, *1/ */
/*   /1*     groupRanks[group], *1/ */
/*   /1*     &physicsGroup); *1/ */
/*   /1* MPI_Comm_create(MPI_COMM_WORLD, physicsGroup, &comm_physics); *1/ */
/*   /1* MPI_Group_rank(physicsGroup, &physicsRank); *1/ */

/*   /1* std::cout << " rank:" << rank *1/ */
/*   /1*   << " group:" << group *1/ */
/*   /1*   << " groupRank:" << physicsRank *1/ */
/*   /1*   << " groupRanks[][]:" << groupRanks[group][rank%physicsNodeSize] *1/ */ 
/*   /1*   << std::endl; *1/ */

/*   /1* int errCode = MPI_Barrier( MPI_COMM_WORLD ); *1/ */
/*   /1* int* driverRanks = new int[nGroups]; *1/ */
/*   /1* for (unsigned int i=0; i<nGroups; i++) *1/ */
/*   /1*   driverRanks[i] = groupRanks[i][0]; *1/ */

/*   /1* MPI_Group_incl( *1/ */
/*   /1*     globalGroup, *1/ */ 
/*   /1*     nGroups, *1/ */
/*   /1*     driverRanks, *1/ */
/*   /1*     &driverGroup); *1/ */
/*   /1* MPI_Comm_create(MPI_COMM_WORLD, driverGroup, &driverComm); *1/ */
/*   /1* MPI_Group_rank(driverGroup, &driverRank); *1/ */
/*   /1* int driverGroupSize; *1/ */

/*   /1* if ( rank%physicsNodeSize == 0) *1/ */
/*   /1* { *1/ */
/*   /1*   MPI_Comm_size(driverComm,&driverGroupSize); *1/ */
/*   /1*   std::cout << " rank:" << rank *1/ */ 
/*   /1*     << " physicsGroup:" << group *1/ */
/*   /1*     << " driverRank:" << driverRank *1/ */
/*   /1*     << " driverRanks[]:" << driverRanks[group] *1/ */ 
/*   /1*     << std::endl; *1/ */
/*   /1*   MPI_Comm_rank(driverComm,&rank); *1/ */
/*   /1* } *1/ */

/*   /1* MPI_Comm_rank(MPI_COMM_WORLD,&rank); *1/ */

/*   /1* LibMeshInit libmesh_init(argc, argv, comm_physics); *1/ */

    

/*   /1* //------- initialize driver *1/ */
/*   /1* libMesh::Parallel::Communicator comm(driverComm); *1/ */
/*   /1* libMesh::Parallel::Communicator physicsComm(comm_physics); *1/ */
/*   /1* AGNOS::Driver agnos( comm, physicsComm, inputfile); *1/ */


/*   //====================================== */
/*   // using split comm instead */

/*   libMesh::Parallel::Communicator comm(MPI_COMM_WORLD); */
/*   libMesh::Parallel::Communicator physicsComm; */

/*   int color = rank / physicsNodeSize ; */
/*   int key = rank % physicsNodeSize; */
/*   std::cout << "rank/nodeSize:" << color << std::endl; */
/*   std::cout << "rank\%nodeSize:" << key  << std::endl; */
/*   comm.split( color, key, physicsComm ); */

/*   std::cout << "commSize:" << comm.size() << std::endl; */
/*   std::cout << "commRank:" << comm.rank() << std::endl; */
/*   std::cout << "physicSize:" << physicsComm.size() << std::endl; */
/*   std::cout << "physicsRank:" << physicsComm.rank() << std::endl; */

/*   LibMeshInit libmesh_init(argc,argv, physicsComm.get() ); */

/*   std::cout << "physicSize:" << physicsComm.size() << std::endl; */
/*   std::cout << "physicsRank:" << physicsComm.rank() << std::endl; */

/*   libMesh::Parallel::Communicator newComm(MPI_COMM_WORLD); */
/*   AGNOS::Driver agnos( newComm, physicsComm , inputfile); */



/*   /1* //------ run driver *1/ */
/*   agnos.run( ); */


/*   /1* MPI_Finalize(); *1/ */
/*   return 0; */
/* } */
