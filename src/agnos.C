
#include "agnosDefines.h"
#include "Driver.h"

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

  //TODO group communicators when requested 
  /* MPI_Comm myComm; */
  /* int mpiSplit = */  
  /*   MPI_Comm_split( MPI_COMM_WORLD, comm.rank(), 0, &myComm); */

  int physicsNodeSize = inputfile("parallel/physicsNodeSize",1);
  std::cout << "physicsNodeSize:" << physicsNodeSize << std::endl;

  std::cout << "comm.size\%physicsNodeSize:" 
    << comm.size()%physicsNodeSize << std::endl;
  if (comm.size()%physicsNodeSize)
  {
    std::cerr << "\n\t ERROR: total number of processors must be a multiple of "
      << "physicsNodeSize.\n" <<std::endl;
    exit(1);
  }


  MPI_Comm physicsComm, driverComm;
  MPI_Group physicsGroup, globalGroup, driverGroup;
  int rank, groupRank; 

  // original group
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm_group(MPI_COMM_WORLD,&globalGroup);

  // physics group
  int group = comm.rank()/physicsNodeSize;
  int nGroups = comm.size()/physicsNodeSize;
  std::cout << "number of groups:" << nGroups << std::endl;

  int** groupRanks = new int*[nGroups];
  for (unsigned int i=0; i<nGroups; i++)
    groupRanks[i] = new int[physicsNodeSize];


  for (unsigned int i=0;i<comm.size();i++)
    groupRanks[i/physicsNodeSize][i%physicsNodeSize] = i;


  std::cout << "-----------------------------\n" ;
  for (unsigned int i=0;i<nGroups;i++)
    for (unsigned int j=0;j<nGroups;j++)
    std::cout << "groupRanks[" << i << "][" << j << "]:" 
      << groupRanks[i][j] <<std::endl;
    std::cout << "-----------------------------\n" ;

  std::cout << "-----------------------------\n" ;
    for (unsigned int j=0;j<nGroups;j++)
    std::cout << "groupRanks[0][" << j << "]:" 
      << groupRanks[0][j] <<std::endl;
    std::cout << "-----------------------------\n" ;

    groupRanks[0][0] = 0;
    groupRanks[0][1] = 1;
    groupRanks[1][0] = 2;
    groupRanks[1][1] = 3;


  MPI_Group_incl(
      globalGroup, 
      group,
      groupRanks[group],
      &physicsGroup);
  int err = MPI_Barrier(MPI_COMM_WORLD);

  MPI_Group_rank(physicsGroup, &groupRank);




  std::cout << "globalRank: " << rank 
    << " group: " << group
    << " groupRank: " << groupRank 
    << std::endl;

  /* LibMeshInit libmesh_init(argc, argv, myComm); */
  
  /* LibMeshInit libmesh_init(argc, argv); */


  //------- initialize driver
  /* AGNOS::Driver agnos( comm, inputfile ); */

  /* //------ run driver */
  /* agnos.run( ); */


  MPI_Finalize();
  return 0;
}


