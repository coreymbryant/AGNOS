#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

// local includes
#include "agnosDefines.h"
#include <mpi.h>

int main( int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  Communicator comm(MPI_COMM_WORLD);

  {
    MPI_Comm subComm;
    int mpiSplit =  
      MPI_Comm_split( MPI_COMM_WORLD, comm.rank(), 0, &subComm);
    LibMeshInit libmesh_init(argc, argv, subComm);

    CppUnit::TextUi::TestRunner runner;
    CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
    runner.addTest( registry.makeTest() );
    runner.run();
  }

  int ierr;
  comm.barrier();
  ierr = MPI_Finalize();
  return 0;
}
