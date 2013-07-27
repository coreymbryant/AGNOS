
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
  libMesh::Parallel::Communicator comm(MPI_COMM_WORLD);


  GetPot inputfile( argv[1] );

  int physicsNodeSize = inputfile("parallel/physicsNodeSize",1);
  /* std::cout << "physicsNodeSize:" << physicsNodeSize << std::endl; */

  if (comm.size()%physicsNodeSize)
  {
    std::cerr << "\n\t ERROR: total number of processors must be a multiple of "
      << "physicsNodeSize.\n" <<std::endl;
    exit(1);
  }


  MPI_Comm physicsComm, driverComm;
  MPI_Group physicsGroup, globalGroup, driverGroup;
  int rank, physicsGroupRank, size; 


  // original group
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

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
  MPI_Comm_create(MPI_COMM_WORLD, physicsGroup, &physicsComm);
  MPI_Group_rank(physicsGroup, &physicsGroupRank);


  std::cout << "globalRank: " << rank 
    << " group: " << group
    << " groupRank: " << physicsGroupRank
    << " groupRanks[][]: " << groupRanks[group][rank%physicsNodeSize] 
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
  MPI_Group_rank(driverGroup, &physicsGroupRank);
  int driverGroupSize;

  if ( rank%physicsNodeSize == 0)
  {
    MPI_Comm_size(driverComm,&driverGroupSize);
    std::cout << "driverRank: " << rank 
      << " physicsGroup: " << group
      << " driverRank: " << physicsGroupRank
      << " driverRanks[]: " << driverRanks[group] 
      << std::endl;
  }


  LibMeshInit libmesh_init(argc, argv, physicsComm);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  /* LibMeshInit libmesh_init(argc, argv); */


  //------- initialize driver
  AGNOS::Driver agnos( comm, driverComm, physicsComm, inputfile);

  /* //------ run driver */
  agnos.run( );


  /* MPI_Finalize(); */
  return 0;
}


