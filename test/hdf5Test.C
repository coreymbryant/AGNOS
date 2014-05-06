

#include "agnosDefines.h"
#include "PhysicsCatenary.h"
#include "PseudoSpectralTensorProduct.h"
// local
#include "ioHandler.h"

// boost
#define BOOST_TEST_MODULE hdf5Test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>


using namespace AGNOS;

BOOST_AUTO_TEST_CASE( ioHandler_user_block )
{
  std::string fileName = "test.h5";

  H5File* h5_file = writeUserBlock( fileName );
  H5File* read_h5_file = readUserBlock( fileName );

  delete h5_file;
  delete read_h5_file;
}

BOOST_AUTO_TEST_CASE( ioHadnler_parameters )
{
  std::string fileName = "test.h5";
  H5File* h5_file = writeUserBlock( fileName );

  // construct parameter object
  unsigned int dimension = 3;
  std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;
  for(unsigned int p=0;p<dimension;p++)
    parameters.push_back(
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter("CONSTANT",1.0,1.0) )
        );


  // write parameters to file
  writeParameters( h5_file, parameters );

  delete h5_file;
  h5_file = readUserBlock( fileName );

  //read in parameters from file
  std::vector<std::shared_ptr<AGNOS::Parameter> > newParameters;
  readParameters( h5_file, newParameters );


  BOOST_REQUIRE( (parameters.size() == newParameters.size() ) );
  for(unsigned int i = 0; i<parameters.size();i++)
  {
    BOOST_REQUIRE( (parameters[i]->type() == newParameters[i]->type()) );
    BOOST_REQUIRE( (parameters[i]->min() == newParameters[i]->min()) );
    BOOST_REQUIRE( (parameters[i]->max() == newParameters[i]->max()) );
  }

  delete h5_file;
}

BOOST_AUTO_TEST_CASE( ioHandler_surrogates )
{

  // order
  std::vector<unsigned int> order;
  order.push_back(1);
  order.push_back(2);
  order.push_back(0);
  order.push_back(0);

  // compute solutions names
  std::set<std::string> computeSolutions ;
  computeSolutions.insert("primal");
  computeSolutions.insert("adjoint");
  computeSolutions.insert("qoi");

  GetPot inputfile = GetPot() ;
  // Set dummy inputs for libmesh initialization
  int ac=1;
  char** av = new char* [ac];
  char* name = new char[27];
  strcpy(name, "test");
  av[0] = name;
  MPI_Init(&ac,&av);
  Communicator comm(MPI_COMM_WORLD);
  LibMeshInit libmesh_init(ac, av, comm.get()) ;
  delete av, name;

  std::shared_ptr<AGNOS::PhysicsModel<T_S,T_P> > physics(
      new AGNOS::PhysicsCatenary<T_S,T_P>(comm,inputfile ) 
      );
  std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;
  for(unsigned int i=0; i<order.size();i++)
    parameters.push_back( 
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter(UNIFORM, 1.0,3.0) )
        );

  std::shared_ptr<AGNOS::SurrogateModel<T_S,T_P> > surrogate (
      new AGNOS::PseudoSpectralTensorProduct<T_S,T_P>(
        comm,
        physics,
        parameters,
        order,
        computeSolutions
        )
      );
  /* surrogate->initialize(); */
  surrogate->build();

  // write surrogate to file
  std::string fileName = "test.h5";
  H5File* h5_file = writeUserBlock( fileName );
  writeSurrogate( h5_file, surrogate );
  delete h5_file;

  //read in surrogate from file
  std::vector<unsigned int> readOrder;
  std::set<std::string> readComputeSolutions ;
  std::vector< std::vector<unsigned int> > readIndexSet ;
  std::map<std::string, LocalMatrix> readCoefficients ;

  h5_file = readUserBlock( fileName );
  readSurrogate( h5_file, readOrder, readComputeSolutions, readIndexSet,
      readCoefficients );
  delete h5_file;

  // check for order
  BOOST_REQUIRE( (order.size() == readOrder.size() ) );
  for(unsigned int o=0;o<order.size();o++)
    BOOST_REQUIRE( (order[o] == readOrder[o]) ) ;

  // check for compute solutions
  BOOST_REQUIRE( (computeSolutions.size() == readComputeSolutions.size() ) );
  std::set<std::string>::iterator sit = computeSolutions.begin();
  for(;sit!=computeSolutions.end();sit++)
    BOOST_REQUIRE( (readComputeSolutions.count(*sit) != 0 ) ) ;

  // check for indexSet 
  std::vector< std::vector<unsigned int> > indexSet = surrogate->indexSet() ;
  BOOST_REQUIRE( (indexSet.size() == readIndexSet.size() ) );
  for(unsigned int i=0; i<indexSet.size(); i++)
  {
    BOOST_REQUIRE( (indexSet[i].size() == readIndexSet[i].size() ) );
    for(unsigned int j=0;j<indexSet[i].size();j++)
      BOOST_REQUIRE( (indexSet[i][j] == readIndexSet[i][j] ) );
  }

  // check for coefficients
  std::map<std::string, LocalMatrix> coefficients 
    = surrogate->getCoefficients();
  BOOST_REQUIRE( (coefficients.size() == readCoefficients.size() ) );
  
  std::map<std::string, LocalMatrix>::iterator cit 
    = coefficients.begin(); 
  for(;cit!=coefficients.end();cit++)
  {
    BOOST_REQUIRE( (cit->second.m() == readCoefficients[cit->first].m()) ) ;
    BOOST_REQUIRE( (cit->second.n() == readCoefficients[cit->first].n()) ) ;
    for(unsigned int i=0;i<cit->second.m();i++)
      for(unsigned int j=0;j<cit->second.n();j++)
        BOOST_REQUIRE( (cit->second(i,j) == readCoefficients[cit->first](i,j) ) );
  }

}
