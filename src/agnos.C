
#include "agnosDefines.h"
#include "Driver.h"


int main(int argc, char* argv[])
{
  if( argc < 2 )
  {
    std::cerr << "\n\t ERROR: must specify AGNOS input file.\n" << std::endl;
    exit(1);
  }

  MPI_Init(&argc,&argv);

  MPI_Comm comm_physics, driverComm;
  MPI_Group physicsGroup, globalGroup, driverGroup;
  int rank, physicsRank, driverRank, size; 

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);



  GetPot inputfile( argv[1] );

  int physicsNodeSize = inputfile("parallel/physicsNodeSize",1);
  /* std::cout << "physicsNodeSize:" << physicsNodeSize << std::endl; */

  if (size%physicsNodeSize)
  {
    std::cerr << "\n\t ERROR: total number of processors must be a multiple of "
      << "physicsNodeSize.\n" <<std::endl;
    exit(1);
  }


  // original group
  MPI_Comm_group(MPI_COMM_WORLD,&globalGroup);

  // physics group
  int group = rank/physicsNodeSize;
  int nGroups = size/physicsNodeSize;

  int** groupRanks = new int*[nGroups];
  for (unsigned int i=0; i<nGroups; i++)
    groupRanks[i] = new int[physicsNodeSize];

  for (unsigned int i=0;i<size;i++)
    groupRanks[i/physicsNodeSize][i%physicsNodeSize] = i;

  MPI_Group_incl(
      globalGroup, 
      physicsNodeSize,
      groupRanks[group],
      &physicsGroup);
  MPI_Comm_create(MPI_COMM_WORLD, physicsGroup, &comm_physics);
  MPI_Group_rank(physicsGroup, &physicsRank);

  std::cout << " rank:" << rank
    << " group:" << group
    << " groupRank:" << physicsRank
    << " groupRanks[][]:" << groupRanks[group][rank%physicsNodeSize] 
    << std::endl;

  int errCode = MPI_Barrier( MPI_COMM_WORLD );
  int* driverRanks = new int[nGroups];
  for (unsigned int i=0; i<nGroups; i++)
    driverRanks[i] = groupRanks[i][0];

  MPI_Group_incl(
      globalGroup, 
      nGroups,
      driverRanks,
      &driverGroup);
  MPI_Comm_create(MPI_COMM_WORLD, driverGroup, &driverComm);
  MPI_Group_rank(driverGroup, &driverRank);
  int driverGroupSize;

  if ( rank%physicsNodeSize == 0)
  {
    MPI_Comm_size(driverComm,&driverGroupSize);
    std::cout << " rank:" << rank 
      << " physicsGroup:" << group
      << " driverRank:" << driverRank
      << " driverRanks[]:" << driverRanks[group] 
      << std::endl;
    MPI_Comm_rank(driverComm,&rank);
  }

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  LibMeshInit libmesh_init(argc, argv, comm_physics);

    

  //------- initialize driver
  libMesh::Parallel::Communicator comm(driverComm);
  libMesh::Parallel::Communicator physicsComm(comm_physics);
  AGNOS::Driver agnos( comm, physicsComm, inputfile);


  //======================================
  // using split comm instead

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



  /* //------ run driver */
  agnos.run( );


  /* MPI_Finalize(); */
  return 0;
}


