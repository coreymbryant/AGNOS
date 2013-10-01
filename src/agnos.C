
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

  /* MPI_Init(&argc,&argv); */

  /* MPI_Comm comm_physics, driverComm; */
  /* MPI_Group physicsGroup, globalGroup, driverGroup; */
  /* int rank, physicsRank, driverRank, size; */ 

  /* MPI_Comm_rank(MPI_COMM_WORLD, &rank); */
  /* MPI_Comm_size(MPI_COMM_WORLD, &size); */

  std::string file_name = argv[1];
  GetPot inputfile( file_name);


  MPI_Init(&argc,&argv);
  Communicator comm(MPI_COMM_WORLD);

  /* int physicsNodeSize = inputfile("parallel/physicsNodeSize",1); */
  /* /1* std::cout << "physicsNodeSize:" << physicsNodeSize << std::endl; *1/ */

  /* if (size%physicsNodeSize) */
  /* { */
  /*   std::cerr << "\n\t ERROR: total number of processors must be a multiple of " */
  /*     << "physicsNodeSize.\n" <<std::endl; */
  /*   exit(1); */
  /* } */


  /* // original group */
  /* MPI_Comm_group(MPI_COMM_WORLD,&globalGroup); */

  /* // physics group */
  /* int group = rank/physicsNodeSize; */
  /* int nGroups = size/physicsNodeSize; */

  /* int** groupRanks = new int*[nGroups]; */
  /* for (unsigned int i=0; i<nGroups; i++) */
  /*   groupRanks[i] = new int[physicsNodeSize]; */

  /* for (unsigned int i=0;i<size;i++) */
  /*   groupRanks[i/physicsNodeSize][i%physicsNodeSize] = i; */

  /* MPI_Group_incl( */
  /*     globalGroup, */ 
  /*     physicsNodeSize, */
  /*     groupRanks[group], */
  /*     &physicsGroup); */
  /* MPI_Comm_create(MPI_COMM_WORLD, physicsGroup, &comm_physics); */
  /* MPI_Group_rank(physicsGroup, &physicsRank); */

  /* std::cout << " rank:" << rank */
  /*   << " group:" << group */
  /*   << " groupRank:" << physicsRank */
  /*   << " groupRanks[][]:" << groupRanks[group][rank%physicsNodeSize] */ 
  /*   << std::endl; */

  /* int errCode = MPI_Barrier( MPI_COMM_WORLD ); */
  /* int* driverRanks = new int[nGroups]; */
  /* for (unsigned int i=0; i<nGroups; i++) */
  /*   driverRanks[i] = groupRanks[i][0]; */

  /* MPI_Group_incl( */
  /*     globalGroup, */ 
  /*     nGroups, */
  /*     driverRanks, */
  /*     &driverGroup); */
  /* MPI_Comm_create(MPI_COMM_WORLD, driverGroup, &driverComm); */
  /* MPI_Group_rank(driverGroup, &driverRank); */
  /* int driverGroupSize; */

  /* if ( rank%physicsNodeSize == 0) */
  /* { */
  /*   MPI_Comm_size(driverComm,&driverGroupSize); */
  /*   std::cout << " rank:" << rank */ 
  /*     << " physicsGroup:" << group */
  /*     << " driverRank:" << driverRank */
  /*     << " driverRanks[]:" << driverRanks[group] */ 
  /*     << std::endl; */
  /*   MPI_Comm_rank(driverComm,&rank); */
  /* } */

  /* MPI_Comm_rank(MPI_COMM_WORLD,&rank); */

  /* LibMeshInit libmesh_init(argc, argv, comm_physics); */

    

  /* //------- initialize driver */
  /* libMesh::Parallel::Communicator comm(driverComm); */
  /* libMesh::Parallel::Communicator physicsComm(comm_physics); */
  /* AGNOS::Driver agnos( comm, physicsComm, inputfile); */


  /* //====================================== */
  /* // using split comm instead */

  /* /1* libMesh::Parallel::Communicator comm(MPI_COMM_WORLD); *1/ */
  /* /1* libMesh::Parallel::Communicator physicsComm; *1/ */

  /* /1* int color = rank / physicsNodeSize ; *1/ */
  /* /1* int key = rank % physicsNodeSize; *1/ */
  /* /1* std::cout << "rank/nodeSize:" << color << std::endl; *1/ */
  /* /1* std::cout << "rank\%nodeSize:" << key  << std::endl; *1/ */
  /* /1* comm.split( color, key, physicsComm ); *1/ */

  /* /1* std::cout << "commSize:" << comm.size() << std::endl; *1/ */
  /* /1* std::cout << "commRank:" << comm.rank() << std::endl; *1/ */
  /* /1* std::cout << "physicSize:" << physicsComm.size() << std::endl; *1/ */
  /* /1* std::cout << "physicsRank:" << physicsComm.rank() << std::endl; *1/ */

  /* /1* LibMeshInit libmesh_init(argc,argv, physicsComm.get() ); *1/ */

  /* /1* std::cout << "physicSize:" << physicsComm.size() << std::endl; *1/ */
  /* /1* std::cout << "physicsRank:" << physicsComm.rank() << std::endl; *1/ */

  /* /1* /2* libMesh::Parallel::Communicator newComm(MPI_COMM_WORLD); *2/ *1/ */
  /* /1* /2* AGNOS::Driver agnos( newComm, physicsComm.get() , inputfile); *2/ *1/ */



  /* /1* //------ run driver *1/ */
  /* agnos.run( ); */
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

