
// local
#include "agnosDefines.h"
#include "PhysicsCatenary.h"
#include "PseudoSpectralTensorProduct.h"
#include "ioHandler.h"

// boost
#define BOOST_TEST_MODULE hdf5Test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace libMesh;
using namespace AGNOS;

BOOST_AUTO_TEST_CASE( ioHandler_constructor )
{
  std::string fileName = "test.h5";

  {
    H5IO h5io( fileName, H5F_ACC_TRUNC );
    BOOST_REQUIRE( (h5io.fileName() == fileName) ) ;
    BOOST_REQUIRE( (h5io.accessFlag() == H5F_ACC_TRUNC) ) ;
    h5io.close();
  }
  {
    H5IO h5io( fileName, H5F_ACC_RDONLY );
    BOOST_REQUIRE( (h5io.accessFlag() == H5F_ACC_RDONLY) ) ;
    h5io.close();
  }

}

BOOST_AUTO_TEST_CASE( ioHandler_user_block )
{
  std::string fileName = "test.h5";

  {
    H5IO h5io( fileName, H5F_ACC_TRUNC );
    h5io.writeUserBlock( );
    h5io.close();
  }
  {
    H5IO h5io( fileName, H5F_ACC_RDONLY );
    h5io.readUserBlock( );
    h5io.close();
  }


}

BOOST_AUTO_TEST_CASE( ioHadnler_parameters )
{
  std::string fileName = "test.h5";

  // construct parameter object
  unsigned int dimension = 3;
  std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;
  for(unsigned int p=0;p<dimension;p++)
    parameters.push_back(
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter("CONSTANT",1.0,1.0) )
        );


  // write parameters to file
  {
    H5IO h5io( fileName, H5F_ACC_TRUNC );
    h5io.writeParameters( h5io.file(), parameters );
    h5io.close();
  }


  //read in parameters from file
  std::vector<std::shared_ptr<AGNOS::Parameter> > newParameters;
  {
    H5IO h5io( fileName, H5F_ACC_RDONLY );
    h5io.readParameters( h5io.file(), newParameters );
    h5io.close();
  }


  BOOST_REQUIRE( (parameters.size() == newParameters.size() ) );
  for(unsigned int i = 0; i<parameters.size();i++)
  {
    BOOST_REQUIRE( (parameters[i]->type() == newParameters[i]->type()) );
    BOOST_REQUIRE( (parameters[i]->min() == newParameters[i]->min()) );
    BOOST_REQUIRE( (parameters[i]->max() == newParameters[i]->max()) );
  }

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

  std::shared_ptr<AGNOS::SurrogateModelBase<T_S,T_P> > surrogate (
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

  std::string fileName = "test.h5";
  // write surrogate to file
  {
    H5IO h5io( fileName, H5F_ACC_TRUNC );
    h5io.writeSurrogate( h5io.file(), surrogate );
    h5io.close();
  }
  std::cout << "post write surrogate" << std::endl;
  //read in surrogate from file
  std::shared_ptr<AGNOS::SurrogateModelBase<T_S,T_P> > readSurrogate;
  {
    H5IO h5io( fileName, H5F_ACC_RDONLY );
    h5io.readSurrogate( h5io.file(), readSurrogate, comm );
    h5io.close();
  }
  std::vector<unsigned int> readOrder = readSurrogate->getExpansionOrder();
  std::set<std::string> readComputeSolutions = readSurrogate->getSolutionNames();
  std::vector< std::vector<unsigned int> > readIndexSet =
    readSurrogate->indexSet();
  std::map<std::string, LocalMatrix> readCoefficients  =
    readSurrogate->getCoefficients();


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

BOOST_AUTO_TEST_CASE( ioHandler_element )
{
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

  // physics
  GetPot inputfile = GetPot() ;
  std::shared_ptr<AGNOS::PhysicsModel<T_S,T_P> > physics(
      new AGNOS::PhysicsCatenary<T_S,T_P>(comm,inputfile ) 
      );


  // order
  std::vector<unsigned int> order;
  order.push_back(1);
  order.push_back(2);
  order.push_back(0);
  order.push_back(0);
  
  // parameters
  std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;
  for(unsigned int i=0; i<order.size();i++)
    parameters.push_back( 
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter(UNIFORM, 1.0,3.0) )
        );

  // compute solutions names
  std::set<std::string> computeSolutions ;
  computeSolutions.insert("primal");
  computeSolutions.insert("adjoint");
  computeSolutions.insert("qoi");

  // surrogate
  std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > > surrogates;
  surrogates.push_back(
        std::shared_ptr<AGNOS::PseudoSpectralTensorProduct<T_S,T_P> >(
          new AGNOS::PseudoSpectralTensorProduct<T_S,T_P>( 
            comm, 
            physics,
            parameters, 
            order,
            computeSolutions )
          ) 
      );
  surrogates[0]->build();


  // construct an element
  std::shared_ptr<AGNOS::Element<T_S,T_P> > elemPtr(
      new Element<T_S,T_P>( parameters,surrogates,physics, 1.0));

  std::string fileName = "test.h5";
  // write element to file
  {
    H5IO h5io( fileName, H5F_ACC_TRUNC );
    h5io.writeElement( h5io.file(), elemPtr);
    h5io.close();
  }

  //read in surrogate from file
  std::shared_ptr<AGNOS::Element<T_S,T_P> > readElemPtr ;
  {
    H5IO h5io( fileName, H5F_ACC_RDONLY );
    h5io.readElement( h5io.file(), readElemPtr, comm, comm );
  }

  // check for weight
  double readWeight = readElemPtr->weight() ;
  BOOST_REQUIRE( (1.0 == readWeight) );

  // check for parameters
  std::vector<std::shared_ptr<AGNOS::Parameter> > readParameters
    = readElemPtr->parameters();
  BOOST_REQUIRE( (readParameters.size() == parameters.size()) ) ;
  
  // check for surrogates
  std::vector<std::shared_ptr<SurrogateModelBase<T_S,T_P> > > readSurrogates 
    = readElemPtr->surrogates();
  BOOST_REQUIRE( (readSurrogates.size() == surrogates.size()) ) ;
  
}

BOOST_AUTO_TEST_CASE( ioHandler_simulation )
{
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

  // physics
  GetPot inputfile = GetPot() ;
  std::shared_ptr<AGNOS::PhysicsModel<T_S,T_P> > physics(
      new AGNOS::PhysicsCatenary<T_S,T_P>(comm,inputfile ) 
      );


  // order
  std::vector<unsigned int> order;
  order.push_back(1);
  order.push_back(2);
  order.push_back(0);
  order.push_back(0);
  
  // parameters
  std::vector<std::shared_ptr<AGNOS::Parameter> > parameters;
  for(unsigned int i=0; i<order.size();i++)
    parameters.push_back( 
        std::shared_ptr<AGNOS::Parameter>(
          new AGNOS::Parameter(UNIFORM, 1.0,3.0) )
        );

  // compute solutions names
  std::set<std::string> computeSolutions ;
  computeSolutions.insert("primal");
  computeSolutions.insert("adjoint");
  computeSolutions.insert("qoi");

  // surrogate
  std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > > surrogates;
  surrogates.push_back(
        std::shared_ptr<AGNOS::PseudoSpectralTensorProduct<T_S,T_P> >(
          new AGNOS::PseudoSpectralTensorProduct<T_S,T_P>( 
            comm, 
            physics,
            parameters, 
            order,
            computeSolutions )
          ) 
      );
  // TODO we have to build this surrogate, otherwise we get memory errors
  surrogates[0]->build();

  // list of elements to mimic activeElem list
  std::list<AGNOS::Element<T_S,T_P> > activeElems;

  // construct an element and split into children
  Element<T_S,T_P> element( parameters,surrogates,physics, 1.0) ;
  std::vector< Element<T_S,T_P> > children = element.split() ;
  std::cout << "children: " << children.size() << std::endl;
  for(unsigned int i=0;i<children.size();i++)
  {
    std::cout << "child: " << i << std::endl;
    std::vector< std::shared_ptr<SurrogateModelBase<T_S,T_P> > > childSurrogates;
    childSurrogates.push_back(
          std::shared_ptr<AGNOS::PseudoSpectralTensorProduct<T_S,T_P> >(
            new AGNOS::PseudoSpectralTensorProduct<T_S,T_P>( 
              comm, 
              children[i].physics(),
              children[i].parameters(), 
              order,
              computeSolutions )
            ) 
        );
    childSurrogates[0]->build();
    children[i].setSurrogates( childSurrogates );
    activeElems.push_front(children[i]) ;
  }

  std::string fileName = "test.h5";
  // write element to file
  {
    H5IO h5io( fileName, H5F_ACC_TRUNC );
    h5io.writeSimulation(activeElems);
    h5io.close();
  }

  //read in surrogate from file
  std::list<AGNOS::Element<T_S,T_P> > readActiveElems;
  {
    H5IO h5io( fileName, H5F_ACC_RDONLY );
    h5io.readSimulation(readActiveElems, comm, comm );
  }


  // check for elements
  BOOST_REQUIRE( (readActiveElems.size() == activeElems.size()) );
  
}
