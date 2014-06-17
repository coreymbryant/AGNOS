#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

// local includes
#include "agnosDefines.h"

using namespace libMesh;

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

#ifndef AGNOS_USING_MPICH
  int ierr;
  comm.barrier();
  MPI_Finalize();
#endif
  return 0;
}
