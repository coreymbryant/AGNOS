
#include <mpi.h>
#include "Driver.h"

using namespace libMesh;

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
  PETSC_COMM_WORLD = MPI_COMM_WORLD ;
  int ierr ;
  ierr = PetscInitialize(&argc, const_cast<char***>(&argv),NULL,NULL);

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int physicsNodeSize = inputfile("parallel/physicsNodeSize",1);
  /* std::cout << "physicsNodeSize:" << physicsNodeSize << std::endl; */

  if (size%physicsNodeSize)
  {
    std::cerr << "\n\t ERROR: total number of processors must be a multiple of "
      << "physicsNodeSize.\n" <<std::endl;
    exit(1);
  }


  /* int globalGroup; */
  // original group
  /* MPI_Comm_group(MPI_COMM_WORLD,&globalGroup); */


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

  /* libMesh::Parallel::Communicator comm(MPI_COMM_WORLD); */
  /* libMesh::Parallel::Communicator physicsComm; */

  /* int color = rank / physicsNodeSize ; */
  /* int key = rank % physicsNodeSize; */
  /* std::cout << "rank/nodeSize:" << color << std::endl; */
  /* std::cout << "rank\%nodeSize:" << key  << std::endl; */
  /* comm.split( color, key, physicsComm ); */

  /* std::cout << "commSize:" << comm.size() << std::endl; */
  /* std::cout << "commRank:" << comm.rank() << std::endl; */
  /* std::cout << "physicSize:" << physicsComm.size() << std::endl; */
  /* std::cout << "physicsRank:" << physicsComm.rank() << std::endl; */

  /* LibMeshInit libmesh_init(argc,argv, physicsComm.get() ); */

  /* std::cout << "physicSize:" << physicsComm.size() << std::endl; */
  /* std::cout << "physicsRank:" << physicsComm.rank() << std::endl; */

  /* /1* libMesh::Parallel::Communicator newComm(MPI_COMM_WORLD); *1/ */
  /* /1* AGNOS::Driver agnos( newComm, physicsComm.get() , inputfile); *1/ */



  /* /1* //------ run driver *1/ */
  /* agnos.run( ); */
  {
    int mpiSplit;
    int physicsColor, physicsKey;

    MPI_Group driverGroup, globalGroup;
    MPI_Comm physicsComm, driverComm;

    // original group
    MPI_Comm_group(MPI_COMM_WORLD,&globalGroup);

    // determine where this proc belongs
    physicsColor = rank / physicsNodeSize ;
    physicsKey = rank % physicsNodeSize;

    // split communicator into physics groups
    mpiSplit =  
      MPI_Comm_split( MPI_COMM_WORLD, physicsColor, physicsKey, &physicsComm);
    LibMeshInit libmesh_init(argc, argv, physicsComm);

    std::cout << "physicsColor:" << physicsColor << std::endl;
    std::cout << "physicsKey:" << physicsKey  << std::endl;


    // create group for driver communicator
    // rank0 from each physics group belongs in this group
    int stride = physicsNodeSize;
    int lastRank = size-stride;
    /* if (lastRank == 0) */
    /*   lastRank = size - 1; */
    int ranges[1][3] = { {0,lastRank,stride} } ;
    std::cout << "ranges: 0," << stride << "," << lastRank
      << std::endl;

    mpiSplit = MPI_Group_range_incl(
        globalGroup,
        1,
        ranges,
        &driverGroup);
    mpiSplit = MPI_Comm_create(MPI_COMM_WORLD, driverGroup, &driverComm);


    AGNOS::Driver agnos( 
        Communicator(driverComm), 
        Communicator(physicsComm), 
        inputfile );

    MPI_Group physicsGroup;;
    MPI_Comm_group(physicsComm,&physicsGroup);
    int groupRank;
    MPI_Group_rank(physicsGroup, &groupRank);
    std::cout << "physicsGroupRank:" << groupRank << std::endl;
    if (groupRank==0)
    {
      MPI_Comm_rank(driverComm,&groupRank);
      std::cout << "driverGroupRank:" << groupRank << std::endl;
    }


    std::cout << "rank:" << rank << std::endl;
    std::cout << "size:" << size  << std::endl;
    agnos.run( );
    
  }

  MPI_Barrier(MPI_COMM_WORLD);

#ifndef AGNOS_USING_MPICH
  ierr;
  comm.barrier();
  MPI_Finalize();
#endif
  return 0;
}

